# 并行、着色与接触求值：方案与评估

本文对应四项工作：宽相位（broad phase）并行化、颜色路径里 `preCompute` 每轮两次的问题、
`stampCached` 的 basis 三明治/下三角重写、以及着色算法的取舍（Welsh-Powell 已实测并回退）。
最后给出对整个并行设施的评估和后续方案。并行设施的改造（`TaskFlow`：Serial /
Independent / Region 三分类 + 共享 worker 预算）已经实现，单独写在
`docs/job-pool-plan.md`。

实测环境：AMD Ryzen 5 2500U（4 物理核 / 8 SMT），15 W，Linux，`-m release`。
测量期间机器上有语言服务器等后台负载（load 1–4），因此**绝对毫秒数会漂移**；
凡是 A/B 结论都取自同一时段内的交替测量，或改用 `perf` 的 instructions / task-clock
作为佐证。基准场景是 `avbd_core_test --bench`：529 个刚体、10 次迭代、约 4400 个接触点。

---

## 0. 结论先行

| 项目 | 状态 | 结果 |
| --- | --- | --- |
| 宽相位并行化 + 真基数排序 | 已实现 | 单独 A/B 串行约 −4%；该阶段不再是最差的一块 |
| `preCompute` 每 manifold 每轮两次 | 已实现（位一致） | 串行 −13%，8 线程 CPU −9%；状态摘要不变 |
| `stampCached` basis 三明治 + 手写下三角 | 已实现 | instructions −11%，串行 wall −8.8%，8 线程 −9.2%，见 §3.4 |
| `preCompute`/`updateDual` 的 C 也进 basis 框架 | 已实现 | instructions −15%，串行 wall −14%，8 线程 −13.5%；`JacobianCache` 每 manifold 少 576 B，见 §3.5 |
| Welsh-Powell 着色 | **已回退 greedy** | 实测颜色数不降反升（Pyramid 4→5），见 §4 |
| 并行设施 | 已改造 | `TaskFlow` 三分类 + 共享 `Scheduler` + region 槽位预算；删除全局 `try_lock`/`wait()`。单 solver 无可测回归，双 solver 吞吐 +24%（4 线程），见 `docs/job-pool-plan.md` |
| 并行设施是否还要重构 | 暂不需要 | 指令数跨线程数完全一致，8 SMT 线程 ~2.3×（三明治后）；瓶颈是 BSP 结构与 SMT/频率。任务图（R2）留到 >8 物理核/多世界 |
| 总体 | — | 529 体基准（当前）：串行 17.19 → 11.50 ms/step（−33%），8 线程 8.48 → 5.77（−32%） |

---

## 1. 基准数据

### 1.1 整体

| threads | 改动前 ms/step | +宽相位/preCompute | +stamp 三明治 | +C basis 框架 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 17.19 | 15.21 | 14.05 | 11.50 |
| 2 | 15.69 | 14.87 | 13.70 | 11.08 |
| 4 | 11.55 | 10.35 | 9.24 | 7.34 |
| 8 | 8.48 | 7.91 | 5.90 | 5.77 |
| 12 | 8.46 | 7.86 | 5.89 | 5.70 |

最后一列是 load≈2 下测的（同一时段的受控 A/B 见 §3.5）；绝对毫秒会随负载漂移。
8 线程约 2.0×（4C/8T + 全核降频；不是代码只用了 2 个核）。
用 `perf stat` 验证过：1 线程与 8 线程的 **instructions 几乎完全相同**
（16.69e9 vs 16.71e9，150 步），说明并行路径没有重复计算；
task-clock 从 13.0 ms/step 涨到 34.9 ms/step，是频率（3.34 → 2.14 GHz）与超线程
IPC（2.50 → 1.49）造成的，不是自旋浪费。

### 1.2 热区（`perf report`，8 线程，C 也进 basis 框架之后）

| 符号 | 占比 |
| --- | ---: |
| `Manifold::stampCached` | 35.8% |
| `Manifold::updateDual` | 15.1%（C 重构前 22.7%） |
| `Solver::iterateBody` 的循环体（Block、力遍历、solve） | 14.6% |
| `Manifold::preCompute` | 12.8%（C 重构前 14.7%） |
| 宽相位（hash + radix + 精确过滤 + pair walk） | ~12% |
| `collideShapes` + `initialize` | ~3% |

现在最大的一块是 `stampCached` 的内存访问和数学本身；接下来性价比最高的是
manifold 布局（SoA / 4 接触点分块）和 `stampCached` 里 S 的跨 body 复用（§6 P4）。

`std::barrier` 的等待只占 0.31% 的采样：在当前规模下，**每颜色一个 barrier 不是主要成本**。
这很关键，因为它意味着“减少颜色数”本身不会带来明显收益；颜色顺序的价值只能体现在收敛质量上。

---

## 2. 宽相位（已实现）

### 2.1 做了什么

1. **每个 body 的准备工作并行**：AABB 半长、层级 `levels[i]`、排序键 `keyed[i]` 都按
   `i` 写入预分配数组（原先 `push_back`，现在是下标写入），不再串行。
2. **修正偷懒**：`Solver::broadPhase()` 原来在栈上建 `detail::JobPool(8)`，完全忽略
   `Solver::threads`；现在用 solver 自己的 pool。
3. **真基数排序**：
   - 候选对打包成 `uint64`（`min<<32 | max`），`radix_std::sort` 才会走基数路径；
     原来 `std::sort` 比较 `std::pair<int,int>`，名字叫 radix 实际是 introsort。
   - 层级哈希键改成 “20 bit × 3 坐标 + 高 4 bit 层级”，一次 8-bit LSD
     基数排序（带 body payload）同时得到层级分段和格内 body 列表。
4. **pair walk 并行**：新增 `JobPool::forBlocks`，把 body 连续分段，每段写自己的
   输出缓冲（容量跨帧保留），最后拼接。排序+去重在后，所以分段只影响中间顺序，
   不影响最终 manifold 创建顺序，也不影响确定性。
5. **阈值**：这些每 body 只有几十 ns 的准备工作用 `kMinBodiesPerCheapDispatch = 512`
   单独设阈值；低于它串行执行，避免为 30 µs 的工作付一次 dispatch。

### 2.2 剩余的串行部分

- `baseCell`（中位数）与 `top` 是两次 O(n) 扫描；
- manifold 的构造/链接必须串行（侵入式链表 + 静态池互斥锁）；
- 精确过滤已经并行。

宽相位实测约 1.7 ms/step（串行），其中构造 manifold、`collide`、`initialize` 不在
`broadPhase` 内。要继续压缩，方向是“批量构造 manifold 后再串行链接”，而不是再加线程。

---

## 3. 颜色路径的 `preCompute`：每轮两次 → 一次（已实现）

### 3.1 问题

每轮每个 manifold 会被两个 body 各访问一次，`Manifold::updatePrimal` 里
`preCompute + stampCached`，于是共享的 C/F、两个 body 的角雅可比都被算了两次。
`updateDual` 还会再算一次 C/F。每轮实际是 **2×preCompute + 1×dual**。

### 3.2 被否掉的方案

- **每轮开头加一个 force 并行 pre-pass**：省一次 preCompute，但所有 stamp 都看
  轮首快照，等于把颜色 GS 降级成 Jacobi，还多一个 barrier。不划算。
- **只算自己那一半雅可比**：省不掉 C/F，收益小、代码分叉多。

### 3.3 采用的方案：dual 预置缓存 + 首次触碰消费

事实链：

1. 上一轮 `updateDual` 结束时算出的 `contactForceCache = penalty_new*C + lambda_new`
   与“下一轮轮首做一次 `preCompute`”逐位相同（同 pose、同表达式、同顺序）。
2. 一个 manifold 的两个 body 颜色必然不同（颜色是独立集），所以每轮**第一个**
   触碰它的 body 一定发生在两个 body 都还没动过的时候；此时 pose 就是轮首 pose。
   中间跑过的其它颜色不会碰这两个 body。
3. 于是第一次触碰可以直接消费 dual 留下的缓存；第二次触碰（更晚的颜色）才重算。
   对静态接触（地面/墙）来说只有一次触碰，于是每轮 0 次 preCompute。

实现只有三处：

```cpp
// Manifold
if (!primalCacheReady) preCompute(alpha);
stampCached(body, block);
primalCacheReady = false;   // 消费掉；下一个颜色会重算
```

`preCompute()` 末尾置 `primalCacheReady = true`，`updateDual()` 末尾也置 `true`
（为下一轮准备），`initialize()` 置 `false`（新接触集，旧缓存作废）。

**安全性**：颜色 pass 之间本来就有 barrier；同一颜色内不会有两个 body 碰同一个
manifold（独立集定义），所以标志的读/写没有竞争。overflow 路径的 phase A
仍然无条件 `preCompute`，语义不变。

**位一致**：默认参数（`jacobianRebuildDistance = 0`）下，`stack_300` 摘要
`0x325f1faf5aab76b6` 改动前后完全相同（当时 86 项断言全过；后续三明治/C 重构
按预期各重锚一次）。
只有开了 `jacobianRebuildDistance > 0` 时，dual 可能复用漂移门限内的世界力臂，
这本来就是那条路径的既定近似。

**实测**（同一时段 A/B，新宽相位 + 当时的着色方案）：

| threads | 改前 ms/step | 改后 ms/step | 变化 |
| ---: | ---: | ---: | ---: |
| 1 | 17.71 | 15.42 | −12.9% |
| 8 | 8.72 | 8.15 | −6.5% wall / −9.4% CPU |

### 3.4 `stampCached`：basis 三明治 + 手写下三角（已实现）

每轮每个 manifold 的两个 body 各调一次 `stampCached`，旧实现每个接触点做三个
通用 3×3 乘法和两次转置。利用 `Jlin = ±basis`：

```
S   = B^T K B                 (线性 Hessian，两侧相同，side-independent)
RS  = R S                     (R = skew(rWorld))
H   = R S R^T                 (角 Hessian，手写 6 个独立分量)
rhsLin = ±B^T F
rhsAng = ±R (B^T F)
```

`K` 是对角阵，所以 `S` 只是三条 basis 行外积的加权和（三明治），9 个元素里只有
6 个独立；`H` 的 6 个分量也用闭环公式手写（已用数值脚本对照 `R S R^T` 验证，误差
5e-15）。对称结果镜像写入，保证 `Block::solveDelta` 的 pivot fallback 仍看到完整矩阵。

受控 A/B（同一对二进制交替运行）：

| 指标 | 旧 stamp | 新 stamp | 变化 |
| --- | ---: | ---: | ---: |
| instructions（150 步） | 16.05e9 | 14.25e9 | **−11.2%** |
| cycles | 6.16e9 | 5.75e9 | −6.6% |
| 1 线程 ms/step（最小） | 14.62 | 13.33 | −8.8% |
| 8 线程 ms/step（最小） | 7.06 | 6.41 | −9.2% |

**注意**：浮点结合顺序变了，`stack_300` 摘要重锚；1/2/3/4/8/12 线程仍逐位一致，
全部断言通过。这和之前“手写 LDLT”那次重锚摘要的性质相同。

### 3.5 `preCompute`/`updateDual` 的 C：进入 basis 框架（已实现）

旧实现为 C 构造两个角雅可比矩阵（每接触点 6 个 cross），再做 4 个矩阵-向量乘；
stamp 改成 basis 框架后这两个矩阵已经没人用了。利用同一组恒等式：

```
C_lin = B (dqA_lin - dqB_lin)
C_ang = B [ (dqA_ang x rA) - (dqB_ang x rB) ]
C     = C0(1-alpha) + B [ (dqA_lin - dqB_lin) + (dqA_ang x rA) - (dqB_ang x rB) ]
```

`preCompute` 和 `updateDual` 共用这个 bracket；`JacobianCache` 只剩 2 个 `float3`
（每 manifold 省 8×72 = 576 B），`driftUsable` 现在只决定是否复用世界力臂，
`jAAng/jBAng` 整个删掉。

受控 A/B（新旧二进制交替，同一时段）：

| 指标 | 旧 C | 新 C | 变化 |
| --- | ---: | ---: | ---: |
| instructions（150 步） | 14.30e9 | 12.14e9 | **−15.1%** |
| cycles | 5.90e9 | 5.16e9 | −12.6% |
| 1 线程 ms/step（最小） | 13.25 | 11.72 | −12% |
| 8 线程 ms/step（最小） | 7.07 | 6.13 | −13.5% |

摘要再次重锚：`stack_300 = 0x31bd6ae6440b4a5a`，
`mixed_240 = 0xf6fe93b6f1d871a7`，`soft-body = 0x262a14e1286d94ea`；
1/2/3/4/8/12 线程仍逐位一致，88 项断言全过。

---

## 4. 着色：Welsh-Powell 已实测并回退 greedy

### 4.1 做了什么、看到了什么

按 Welsh-Powell（度数降序、index 升序打破平局、逐颜色建类）实现并跑了全部场景。
颜色数（120 步稳定后）：

| 场景 | greedy | Welsh-Powell | overflow |
| --- | ---: | ---: | ---: |
| Pyramid（136） | **4** | 5 | 0 |
| Pyramid XL（1653） | 6 | 6 | 0 |
| Stack（10） | 2 | 2 | 0 |
| Rope（19） | 2 | 2 | 0 |
| Bridge（88） | 5 | 5 | 0 |
| Soft Body（192） | 8 | 8 | 6 |

Welsh-Powell 没有减少颜色数，Pyramid 还多一种（多一个 barrier/轮）。原因是
Pyramid 的建体顺序恰好接近该图的 4 色解，按度数重排反而打乱它；Soft Body 的
8 色 + overflow 是图本身稠密。`std::barrier` 等待只占 0.31% 采样，所以“颜色数”
本身也不是主要成本——颜色顺序的价值只可能在 GS 收敛质量上，而那需要按
“稳定后穿透 / 达到同一阈值所需轮数”来评，不是按颜色数。

### 4.2 结论

已**回退到原来的 body-order greedy**，保持现有行为和摘要。后续若还要动着色：

1. 接触求值的 P0/P1 已完成；下一步是 §6 P4 的 manifold 布局；
2. 如果要做收敛 A/B，用 Pyramid / Pyramid XL / Soft Body / Bridge 量穿透和轮数；
3. 真要更少颜色，用 **DSATUR**（动态选饱和度最大的点），而不是 Welsh-Powell；
   但按 0.31% 的 barrier 占比，预期总收益 <1%。

---

## 5. 并行设施评估（是否有局限？）

### 5.1 现状清单（`TaskFlow` 改造后）

- `TaskFlow`（`taskflow.hpp`）：`Serial` / `Independent` / `Region` 三类阶段，线性执行。
- `JobPool::forCount` / `forBlocks`：调用线程跑 chunk 0，其余 detach，用 `std::latch`
  做**每次分派**的 join；不再有 `pool.wait()`。
- `JobPool::runRegion` + `LoopSync`：一个 step 的迭代集只 dispatch 一次；
  每个参与者跑完整循环，颜色/对偶阶段之间用 `std::barrier` 同步。
- 进程级 `Scheduler` 共享 `BS::thread_pool`；region 用 `SlotBudget` 预留
  `participants-1` 个 worker，两个 solver 可以按预算重叠，抢不到预算时等待，
  不再静默串行，也没有全局锁。
- 每阶段的 break-even 阈值（48 body / 256 force）；region 还有 `inlineWhen`
  的小场景内联。

### 5.2 实测判断

- **没有重复计算**：instructions 跨线程数不变。
- **barrier 不贵**：0.31% 采样。
- **dispatch 有固定成本**：每步约 5–6 次 fork-join，按代码注释的 45 µs/次估算，
  约 0.2–0.3 ms/step；在 8 线程的 8 ms 里占 3% 左右。
- **8 线程约 2.3×**（三明治后）是 4C/8T + 全核降频的结果；4 线程约 1.8× 也符合该 CPU 的功耗墙。

### 5.3 结构性局限（按影响排序）

1. **每个 step 仍是多次 fork-join**（改造成 `TaskFlow` 后结构更清楚，但阶段数没变）：
   宽相位内部 2 次、精确过滤、warmstart forces、warmstart bodies、finish bodies。
   把它们搬进一个持久并行区是 P2，预计约 3%。
2. ~~**`runLoop` 是全局独占的**~~：已由 `TaskFlow` + `SlotBudget` 解决。第二个
   solver 不再静默串行；region 按预算重叠。`sharedPool().wait()` 已删除（它与阻塞
   在 barrier 上的参与者不能共存）。细节和实测见 `docs/job-pool-plan.md`。
3. **BSP + 每颜色 barrier**：所有线程在每个颜色末尾对齐，一个颜色里最慢的 body
   决定相位长度。工作量按“body 数”而不是“接触数/耗时”均衡，静态大块、
   多接触的 body 会拖尾。颜色数 ≤ 8 是对 barrier 数量的硬上限，也导致 overflow
   的 3 个 barrier + 全 force 扫描。
4. **cursor 分配器**：`LoopSync::forItems` 用共享原子游标，粒度 `count/(4*参与者)`；
   对很小的颜色（比如 20 个 body）会退化成串行，且每相位重置游标。
5. **确定性约束**：不允许跨线程浮点归约，所以自适应穿透统计、迭代数决策只能串行；
   这些都很便宜，可以接受。
6. **内存布局**：`Manifold` 约 1.4 KB AoS，`stampCached` 每个接触点做 3 个 3×3
   乘法和若干转置；现在它是最大单项（36.7%）。这不是并行设施的问题，
   但它是当前性价比最高的优化点。

### 5.4 结论

当前规模下**不需要为了并行而重构**：设施能正确、确定地跑，扩展上限主要来自
4C/8T SMT 和 BSP 负载不均，而不是全局池本身。只有当目标是
“>6–8 物理核”或“多 solver / 后台步进”时，才值得上 §6 的 R1/R2。
短期收益最大的顺序是 P1（数学）→ P2（单并行区）→ P3（着色收敛验证）。

---

## 6. 另一个方案（如果要继续重构）

### R0（推荐）：增量演进，不动架构

- ~~**P1 `stampCached` 的数学**~~：已完成（§3.4）。
- ~~**P0 把 C 也改成 basis 框架**~~：已完成（§3.5），instructions −15%，
  每 manifold 省 576 B。
- **P2 一个 step 一个并行区**：`TaskFlow` 已经把 step 写成线性阶段序列，下一步是
  把 broadphase / warmstart / finish 的 Independent 阶段也放进一个持久 region，
  用 `sync.forItems` 分派，省掉每步 4–5 次 dispatch（约 3%）。注意 region 内不能再
  调用阻塞式 `forCount`，所以宽相位的 chunk 缓冲要改成按 worker 索引。
- **P3 着色收敛 A/B**：用 Pyramid / Pyramid XL / Soft Body / Bridge 量
  “达到同一穿透阈值所需轮数”，再决定 Welsh-Powell / DSATUR / greedy。
  顺带把 overflow 的 phase A 从“扫全 force”改成“只扫 overflow 相关 manifold 列表”。
- **P4 manifold 布局**：`Contact`/`JacobianCache`/`contactForceCache` 拆成 SoA
  或按 4 接触点分块，减少 1.4 KB 结构在 preCompute/dual 两次遍历时的 cache miss。

### R1：按岛（island）分解 + 岛内颜色

- 复用 `islandSleep` 的 union-find 把约束图分成连通分量；岛之间零依赖。
- 小岛直接整岛交给一个 worker 顺序 GS（没有全局 barrier，locality 更好）；
  大岛（例如整面墙/大堆）内部仍用颜色并行，避免单岛串行。
- 调度：所有岛作为一个 work queue，大岛拆成颜色任务，小岛合并成批；
  需要一个“大岛阈值”和确定性的任务编号。
- 收益：去掉“所有线程等最慢颜色”的全局同步，特别适合多个独立堆叠的场景；
  代价：实现复杂，且单个巨大岛仍受颜色 barrier 限制。
- 确定性：每个岛内部保持固定的颜色/下标顺序，岛间无共享写，结果与全局颜色
  顺序在该岛上的子序列一致。

### R2：依赖计数（task graph）代替颜色 barrier

- 每个 body 一个任务，前驱是“所有约束邻居中排在前面的 body”；前驱计数归零即可运行。
- 用一个全局原子计数器 + 就绪队列（persistent workers），没有阶段 barrier。
- 颜色只作为确定性的 tie-break（保证同一次运行中执行顺序固定）。
- 收益：消除 barrier 和颜色阶段，负载自适应；代价：调度开销、调试难度、
  以及必须严格证明“就绪判定”与 GS 顺序等价才能保持位一致。

---

## 7. 验证清单

- [x] `xmake run avbd_core_test`：88 checks / 0 failures（含 1/2/3/4/8/12 线程摘要一致，
      以及两个 solver 并发步进的摘要一致）
- [x] `stack_300` 摘要：preCompute 优化保持 `0x325f1faf5aab76b6` 位一致；
      三明治重写后、C 重构后各重锚一次（当前 `0x31bd6ae6440b4a5a`，预期内）
- [x] `avbd` GDExtension 目标编译通过
- [ ] 着色算法（greedy/WP/DSATUR）的穿透与收敛轮数 A/B（已回退 greedy，暂不做）
- [ ] ASan / TSan 各跑一遍（`primalCacheReady` 的跨阶段可见性、`TaskFlow` 的
      latch/barrier/condvar 都值得用 TSan 确认一次）
- [ ] Pyramid XL / Soft Body / Bridge 的 `--bench-sweep` 在空闲机器上复测

---

## 8. 备注：测量方法

机器上长期有语言服务器（~100% CPU）和浏览器，绝对时间不可靠。本文的 A/B 采用：

- 交替重建、同一时段连续跑，取多次最小值；
- `perf stat -e instructions` 判断“工作是否重复”；
- `perf report` 看占比变化（占比比绝对毫秒稳）；
- 需要精确阶段成本时临时在相位内加 `steady_clock` 累加打印，测完即删。

后续做性能结论时，建议固定 `taskset`、关闭后台索引，并把
`--bench` 的 warmup/timed 步数写进脚本，避免把噪声当回归。

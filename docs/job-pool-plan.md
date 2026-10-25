# TaskFlow：按同步类型分类的执行设施（已实现与复评）

本文记录 `runLoop` 全局 `try_lock` + `sharedPool().wait()` 的替代实现，以及在实现
之后对方案的重新评估。背景与各阶段性能见 `docs/parallel-primal-plan.md`。

---

## 1. 旧实现的问题（已验证）

- `runLoop` 用 `try_lock(loopMutex())` 独占整个进程的持久循环；抢不到锁的 solver
  **静默退化成串行**，没有任何信号。
- `runLoop` / `forCount` / `forBlocks` 都用 `sharedPool().wait()` 做 join，这是
  “等全局队列清空”，不是“等本次分派”。两个 solver 会互相等；更严重的是：
  **`wait()` 与阻塞在 barrier 上的 region 参与者不能共存**——`wait()` 等
  `tasks_running == 0`，而参与者正卡在 barrier 上，这就是潜在死锁。
- `forCount` 不受 `loopMutex` 保护，“独占”只对别的 `runLoop` 成立。
- `LoopSync` 的 barrier 假设恰好 `workers` 个参与者常驻，这个假设依赖共享队列里
  没有别的 solver 的任务，代码没有强制。

---

## 2. 已实现的系统

### 2.1 三个任务类别（`src/avbd/taskflow.hpp`）

```cpp
enum class TaskKind { Serial, Independent, Region };
```

| 类别 | 语义 | 执行 |
| --- | --- | --- |
| `Serial` | 必须在调用线程单独运行（链表修改、删除、着色） | 直接调用 |
| `Independent` | `fn(i)`，item 之间无同步，谁算哪项都不影响结果 | `JobPool::forCount`，每次分派一个 latch join |
| `Region` | `body(worker, sync)`，阶段之间用 `sync.arriveAndWait()` 同步 | `JobPool::runRegion`：预留参与者预算、detach helpers、latch join |

`TaskFlow` 是**线性**阶段序列（不是 DAG）：物理步本身就是线性的，分类的价值在于
每个阶段被允许怎样使用 worker，而不是重排。`independent` 的 `count()` 在
**运行时**求值（更新集/force 列表在同一 step 的前面阶段才确定）；`region` 还接受一个
`inlineWhen()` 谓词，小场景直接走 `NullSync` 内联，避免为几十个 body 唤醒 worker。

### 2.2 调度器与 join（`src/avbd/job_pool.hpp`）

- `Scheduler::shared()`：进程级 `BS::thread_pool`，惰性创建；串行 solver 不创建线程。
- `SlotBudget`：一个 all-or-nothing 的参与者预算（mutex + condvar）。
  一个 region 需要 `helpers = participants - 1` 个 pool worker，**先一次性拿满**
  再 detach；拿不满就等。没有“持有部分槽位再去等其余槽位”，因此两个 region
  不会互相死锁。
- `runChunks`：`forCount`/`forBlocks` 的统一实现。调用线程跑 chunk 0，其余 chunk
  detach 到 pool，用 `std::latch` join；**不再使用 `pool.wait()`**。
- `runRegion`：`participants = threadCount()`，`helpers = participants - 1`；
  `SlotBudget::acquire(helpers)` → detach `helpers` 个参与者 → 调用线程跑
  participant 0 → `latch.wait()` → release。`runRegionInline` 用 `NullSync`。

为什么不会“没有空闲 worker 导致死锁”：

- region 只 detach `helpers = participants - 1 ≤ capacity - 1` 个任务，pool 至少有
  一个 worker 不在 region 里；调用线程本身是 participant 0，始终在推进。
- 非 region 任务（`forCount`/`forBlocks`）不阻塞，最终会完成；它们可能让 region
  的某个 helper 晚一点被调度，但不会让 region 永久停住。
- 两个 region 的 `helpers` 之和超过容量时，后到的 region 在 `SlotBudget` 上等待，
  而不是静默串行。

### 2.3 Solver 的 step 变成 flow

`Solver::step()` 现在显式构造：

```cpp
serial      broadPhase();                 // 内部仍是 independent 阶段
serial      forceCount = warmstartForces();
serial      selectUpdateBodies();
independent warmstartBody(i)              // count 运行时取 warmstartOrder.size()
serial      computeIterations(target, newton);
region      iterateBody(...)              // 小场景 inlineWhen -> NullSync
independent finishBody(i);
```

`iterate` 被拆成 `iterateBody<Sync>`（模板，`LoopSync` / `NullSync` 共用一份循环体）；
`warmstartBodies` / `finishVelocities` 拆成逐 body 的 `warmstartBody(i)` /
`finishBody(i)`，由 flow 的 Independent 阶段调度。`converged` 原子和
`prevResidual` 由 flow 的 region lambda 捕获，所有参与者共享，语义与旧的
`runLoop` 内共享变量一致。

---

## 3. 正确性与确定性

- 改造当时 `stack_300` 摘要 `0xcd475375c89b961b`、`mixed_240` 摘要
  `0x22810f20f393dc98` 与改造前（greedy + 三明治）完全一致；之后的 C basis
  重构按预期又重锚了一次（当前见主文档 §3.5），与本设施无关。
- 1/2/3/4/8/12 线程摘要仍全部一致。
- 新增 `test_concurrent_solvers`：两个 solver 各 4 线程、从两个 host 线程同时跑
  240 步，两个摘要都等于串行参考。这个用例在旧设计下不会死锁（会静默串行），
  但会暴露 `wait()` 的全局 drain 行为；现在它同时覆盖并发与确定性。
- 代码里保留两条规则并写入注释：
  1. region 内不得再 `forCount`/`forBlocks`（参与者阻塞等待新任务）；
  2. region 内不得再起 region（参与者在 barrier 上）。
  两者在当前 step 结构里都不会出现。

---

## 4. 实测

机器：Ryzen 5 2500U（4C/8T），有后台负载，取同一时段交替运行的最小值。

### 4.1 单 solver：没有可测回归

| threads | 旧设施 ms/step | TaskFlow ms/step |
| ---: | ---: | ---: |
| 1 | ~15.9 | ~16.0 |
| 8 | ~8.85–9.2 | ~8.9–10.0 |

流本身的构造（约 7 个 `std::function` 阶段）在 step 的毫秒量级里不可测；
`forCount` 的 chunk 数从“每 worker 1 块”变成“每 worker 4 块”，反而改善了一点
负载均衡。

### 4.2 两个 solver：旧设施串行、新设施按预算重叠

`--bench-two`：两个 529 体基准场景，各跑 120 步，两个 host 线程同时驱动。

| threads（每个 solver） | 旧设施 ms/step/solver | TaskFlow ms/step/solver | 说明 |
| ---: | ---: | ---: | --- |
| 4 | 19.43 | **15.63** | 两个 region 各 3 helpers，6/8，可重叠 |
| 8 | 17.24 | **15.81** | 各 7 helpers，14/8，region 串行，仅非 region 阶段重叠 |
| 1 | — | 20.58 | 两个 host 线程各自串行，受 4 物理核限制 |

4 线程时双 solver 吞吐提升约 **24%**（19.43 → 15.63），8 线程时约 **9%**；
这正是“按预算重叠”而不是“静默串行”的差别。单 solver 的 `--bench` 基线（4 线程
9.98、8 线程 8.80）用于对照，说明重叠不是靠牺牲单 solver 性能换来的。

---

## 5. 重新考虑后的结论

### 5.1 这个系统证明了什么

- 物理步可以干净地分成 **Serial / Independent / Region** 三类；只有 Region 需要
  “所有参与者到齐”的同步，Independent 阶段天然可跨 solver 重叠。
- 旧 `try_lock` 想解决的“队列里有别人的阻塞参与者”问题，用 **按分派 latch +
  region 槽位预算** 解决更准确：既没有静默串行，也没有全局 drain。
- 分类本身不会减少 Region 内部的工作或 barrier；它的价值是**多 solver 并发**和
  **可组合**（以后可以在 Independent 阶段之间做融合，或把 Region 换成任务图）。

### 5.2 现在还不值得上完整 task graph（原方案 C）

- 当前 step 是线性的；Region 是一个不可拆的 barrier 序列。完整任务图的主要收益
  是去掉 region barrier 和按依赖计数调度，但它需要：
  - 每轮/每颜色的依赖计数与就绪队列；
  - 把 per-worker 的 Chebyshev 状态改成 per-round / per-body；
  - 细粒度任务的调度开销（529 体、4C/8T 上很可能净亏）。
- 这些在当前规模没有回报。`TaskFlow` 已经把“哪些工作可以并行、哪些必须 fence”
  这件事抽象出来了；未来要换调度器，只需替换 `JobPool::runRegion` 内部，Solver
  代码不用动。

### 5.3 仍然存在的限制（按影响）

1. **每个 step 仍有多次 fork-join**（宽相位内部 2 次、精确过滤、warmstart forces、
   warmstart bodies、finish bodies）。把它们搬进一个持久并行区是下一步的 P2，
   收益约 3%（见主文档 §6）。
2. **Region 内部是 BSP**：每颜色一个 barrier，工作量按 body 数而非接触数均衡。
   这是算法层的限制，不是设施层的。
3. **`SlotBudget` 非 FIFO**：`notify_all` + 数量谓词，连续的小 region 可能让大
   region 等更久；当前 region 少且短，不构成问题。需要时加 ticket queue。
4. **host 调用线程不计入预算**：多个 solver 并发时，region 参与者 = solver 数 ×
   (threads-1) 受预算约束，但各 solver 的调用线程额外存在；这是有界的短暂超订。
5. **禁止嵌套**：region 内不能再分派阻塞任务。用 debug 断言可以兜住
   （`BS::this_thread::get_pool()` 检查），当前没加。

### 5.4 建议

- **保留** 现在的 `TaskFlow` + `Scheduler` + `SlotBudget`。它直接删掉了全局锁和
  全局 `wait()`，并让多 solver 从“静默串行”变成“按预算重叠”。
- 主文档 §6 P0（C 也改成 basis 框架）已完成：instructions −15%，单/双 solver
  都没有回归；调度器不是当前单 solver 的瓶颈。
- 只有当目标变成“>8 物理核 / 后台物理 / 大量并发世界”时，才考虑把 Region 内部
  换成依赖计数任务图；那时 `TaskFlow` 的类别划分正好是迁移边界。

---

## 6. 验证

- `xmake run avbd_core_test`：**88 checks / 0 failures**（新增 2 项并发 solver 摘要）。
- `--bench-two` 可用；两个 solver 的摘要与串行一致。
- 旧/新设施的 A/B：同一核心测试源码分别链接旧 `job_pool`/`solver` 与新实现；
  数字见 §4。
- 待办：TSan 跑一遍 `test_concurrent_solvers`（latch/barrier/condvar 的可见性）；
  ASan；在空闲机器上复测 §4。

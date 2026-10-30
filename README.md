# Godot-AVBD

![Test](https://github.com/Casinius/godot-avbd/actions/workflows/c-cpp.yml/badge.svg)


## What is AVBD?

AVBD stands for Augmented Vertex Block Descent. It is a position-based physics framework built on top of VBD, using the augmented Lagrangian method to scale the constraint force, which is about as basic as it gets in position-based physics.

It can handle a lot of rigid bodies, because each body only ever solves a small Hessian block of its own, and because a contact between a huge body and a tiny one is just another constraint to it, nothing special.

## What is the benefit?

It does have some. First, it can handle more objects, and more accurately, than physics engines built on the impulse projection method.

And in the future it is entirely possible to move all the computation to the GPU; companies like NVIDIA have already done exactly that with methods of this family.

## How does AVBD work?

### Core formula

AVBD is based on a Hessian-based Newton integration method (extended, of course). From it we get: the force formula

$$
H_i \, \Delta x_i = f_i
$$

and the Hessian iteration formula

$$
H_i = \frac{M_i}{\Delta t^2} + \sum_j H_{ij}
$$

where $M_i$ is the body's mass matrix, and the $H_{ij}$ are the terms every constraint touching body $i$ stamps in. Transform the first formula into this:

$$
\Delta x_i = H_i^{-1} f_i
$$

Each body solves only its own little system. That is the "block" in Vertex Block Descent.

### Physics basic

As we know, real-time simulation pursues minimum energy. So we have this:

$$
E_j(x) = \frac{1}{2} k_j \left( C_j(x) \right)^2
$$

$k_j$ is some kind of stiffness representative, and $C_j$ is one of those "penalty functions", which can tell you how close you are to the physical constraint. Zero means satisfied.

### AVBD magic!

So usually here, XPBD uses simple Lagrange multipliers to get things in position. AVBD uses Lagrange multipliers too, but the augmented kind, the one that can achieve almost infinite stiffness.

So we have:

$$
E_j(x) = \frac{1}{2} k_j \left( C_j(x) \right)^2 + \lambda_j C_j(x)
$$

$\lambda_j$ is the Lagrange multiplier! And here is the thing worth understanding: in a pure penalty method, the only way to enforce a hard constraint is to crank $k_j$ towards infinity. But the constraint oscillates at a frequency of about $\sqrt{k_j / m}$, and the moment that gets close to the step rate, the solver stops converging: the error grows, so you raise the penalty, so the error grows faster. A trap. The multiplier is the way out — the constraint force is carried by $\lambda_j$, not by $k_j$ alone, so $k_j$ never has to reach crazy values.

### The loop

One time-step runs like this: broad phase finds the contacts, every body gets warm-started, then the solver iterates, and each iteration has two passes.

The primal pass: every body assembles its own little system $H_i \, \Delta x_i = f_i$, in which every constraint touching it stamps in its force $F_j = k_j C_j(x) + \lambda_j$, and solves it for a position update. Bodies that share no constraint never read or write each other, so they can all solve at once — this is where graph colouring comes in.

The dual pass: every constraint, with all body poses frozen, updates its own state:

$$
\lambda_j \leftarrow k_j C_j(x) + \lambda_j
$$

so the multiplier accumulates the constraint force iteration by iteration. The penalty itself adapts too:

$$
k_j \leftarrow \min\left( k_j + \beta \left| C_j(x) \right|, \; k_j^{*} \right)
$$

a violated constraint grows its penalty, a satisfied one leaves it alone, and $k_j^{*}$ is a cap (a real one — past it, no step size survives). And at the start of a step the multipliers from the last step are reused, decayed a little, $\lambda_j \leftarrow \alpha \gamma \lambda_j$: warm starting, so a stack that was already still does not rebuild its contact forces from zero every frame.

Contacts get two extras worth naming: the normal force may only push, never pull, and friction lives on a cone, $\left\| F_{\mathrm{t}} \right\| \le \mu \left| F_{\mathrm{n}} \right|$, clamped back into the cone every dual pass.

### So what's next?

That is really the whole story: iterate the two passes until the energy settles, derive the new velocities from the positions, next frame. This repo runs the whole thing multithreaded and proves with bit-exact digests that the threaded result equals the serial one.

The broad phase is a **Morton LBVH with refit** (`src/avbd/bvh_broadphase.*`), ported back from `main`. The repo used to claim the hierarchical hashed grid that replaced it was faster; that comparison crossed branches that also differed by the Eigen-to-linalg port and a full solver rewrite, so it did not hold. Isolated - same swept AABBs, same downstream work, bit-identical results - the LBVH is 1.5-6x faster in the broad phase (3.366 -> 0.953 ms on 1654 bodies) because the hash walks 27 cells at *every* level up to the coarsest, so one large ground body inflates the candidate list 6-9x. See [`tools/baseline/README.md`](tools/baseline/README.md) for what is actually measured (per-phase step cost, the penetration distribution on a 420-box pyramid, raw suite output) and the known limitations.

This is not one of those tutorial articles, so that is all I want to show you. If you are a physics engine dev, you know what I know.

[If you really really want to know](v2-6e8c45cac4eb0cacfc6b952955f36cae_1440w.jpg)

---

# Godot-AVBD（中文版）

## AVBD 是什么？

AVBD 全称 Augmented Vertex Block Descent（增广顶点块下降），是一套基于位置的物理框架：底层是 VBD，在其上引入增广拉格朗日（Augmented Lagrangian）方法来放大约束力。这套思路在 position-based 物理里属于最基础的一类。

它能撑住大规模刚体仿真，原因有二：每个刚体实际求解的只是自身对应的一小块海森矩阵；同时，大刚体与小刚体之间的接触对它而言也只是一条普通约束，不需要任何特殊处理。

## 有什么好处为啥要换他

好处是实打实的。首先，与基于冲量投影法的物理引擎相比，它能处理更多物体，精度也更高 , 其次他能处理超大刚度比和质量比的刚体约束 ,静稳定性更高 ,对墙壁等的组合物体更好

其次，这套方法天然适合向 GPU 迁移——把计算全部搬到 GPU 上完全可行，NVIDIA 已经用同一族的方法做到了。

在我的9750H上 完全清空负载 然后运行4个pyramid交替掉落 , 最低帧数可以保证24到30fps , 实际平均帧率会更高 , 受限于跨平台确定性无法使用simd , 绑核等操作 , 这个性能是还算过得去的

## AVBD 是怎么工作的？

### 核心公式

AVBD 建立在分析力学的牛顿积分法之上（当然做了扩展）,用黑森矩阵作为核心。由此得到两个公式：受力公式

$$
H_i \, \Delta x_i = f_i
$$

以及海森迭代公式

$$
H_i = \frac{M_i}{\Delta t^2} + \sum_j H_{ij}
$$

其中 $M_i$ 是刚体的质量矩阵，$H_{ij}$ 则是每条与该刚体相关的约束所贡献的项。对第一个公式做变形：

$$
\Delta x_i = H_i^{-1} f_i
$$

每个刚体只求解自己的这个小系统——"顶点块下降"里的"块"指的就是这个。

### 物理基础

我们知道，实时物理模拟的本质是最小化能量，于是有：

$$
E_j(x) = \frac{1}{2} k_j \left( C_j(x) \right)^2
$$

$k_j$ 是刚度系数（罚参数），$C_j$ 就是所谓的"罚函数"，衡量当前状态离约束满足有多近，取零表示约束已满足。

### AVBD 的魔法

通常在这一步，XPBD 会引入简单的拉格朗日乘子来处理约束。AVBD 同样使用乘子，但用的是增广版本——这正是它能逼近无限刚度的原因。

于是能量泛函变成：

$$
E_j(x) = \frac{1}{2} k_j \left( C_j(x) \right)^2 + \lambda_j C_j(x)
$$

$\lambda_j$ 就是拉格朗日乘子。这里有一个关键问题值得展开：纯罚函数方法里，想要硬约束，唯一的手段是把 $k_j$ 推向无穷大。但约束会以大约 $\sqrt{k_j / m}$ 的频率振荡，一旦这个频率接近仿真步频，求解器直接不收敛——误差变大，于是加大罚参数，误差变得更大，形成恶性循环。乘子正是打破循环的手段：约束力由 $\lambda_j$ 承担，而不是全部压在 $k_j$ 上，$k_j$ 也就无需取到病态量级。

### 循环

一个时间步的流程：宽相检测找出接触点，对各刚体做热启动，然后求解器开始迭代，每轮迭代分两个 pass。

primal pass：每个刚体组装自己的小系统 $H_i \, \Delta x_i = f_i$，每条相关约束把自身约束力 $F_j = k_j C_j(x) + \lambda_j$ 叠加进去，解出位置增量。不共享约束的刚体之间没有数据依赖，可以完全并行——图染色就是为此服务的。

dual pass：在所有刚体位姿冻结的前提下，每条约束独立更新自身状态：

$$
\lambda_j \leftarrow k_j C_j(x) + \lambda_j
$$

乘子就这样逐轮累积约束力。罚参数同时自适应调整：

$$
k_j \leftarrow \min\left( k_j + \beta \left| C_j(x) \right|, \; k_j^{*} \right)
$$

未被满足的约束增大罚参数，已满足的保持不变；$k_j^{*}$ 是硬上限——超过它，任何步长都无法保证收敛。另外，每个时间步开始时会复用上一步的乘子，仅做少量衰减：$\lambda_j \leftarrow \alpha \gamma \lambda_j$。这就是热启动（warm start）：已经稳定的堆叠不必每帧从零重建接触力。

接触约束还有两个细节值得单独说明：法向力只能推不能拉；摩擦约束位于锥面上，$\left\| F_{\mathrm{t}} \right\| \le \mu \left| F_{\mathrm{n}} \right|$，每次 dual pass 都会把它重新夹回锥内。

### 写在最后

到这里整个算法已经完整：两个 pass 反复迭代直到能量收敛，再由位置差分出新的速度，进入下一帧。这个仓库把整套流程跑在多线程上，并用位级一致的 digest 验证了多线程结果与串行实现完全相同。

宽相默认是 **分层哈希网格**（`src/avbd/spatial_hash.*`），`broadPhaseMode = 1`（`physics/avbd/broad_phase_lbvh=true`）可切到 **Morton LBVH + refit**（`src/avbd/bvh_broadphase.*`）。两个实现吃同一份 swept AABB，产出的线序也相同（最终 digest 逐位一致：`stack_300 = 0xfbac4f08cb1c977c` 等）。

哈希之前的问题不是“树更好”，而是候选对放大：每个刚体把所有 27 个邻格都发出去、同层配对从两侧各发一次，把精确过滤和排序的活全留给了下游——Pyramid XL 上 45603 个候选对，而真实 AABB 重叠只有 4868 个。现在的 walk 做了三件事：

* 同层配对只由较小 body index 发一次；
* 在 walk 内先做 swept AABB 精确重叠测试，再写入候选（下游的 layer / constraint 过滤照旧）；
* 空层直接跳过（一个大地面会把 `top` 抬高，但不会填充中间层）。

结果是候选对与 LBVH 完全相同，Pyramid XL 的整步从优化前 **10.019 ms → 8.70 ms**（`--bench-scene "Pyramid XL" --threads 8`，120 步平均；同一台机器上 LBVH 约 8.35 ms，落在噪声内）。隔离看宽相索引本身，哈希仍比 LBVH 慢（1.18 vs 0.59 ms），但下游已经不再为 9 倍假阳性买单：

| 场景 | 刚体 | 候选对 hash / lbvh | 宽相 hash (其中索引) | 宽相 lbvh (其中索引) | 整步 hash → lbvh |
|---|---|---|---|---|---|
| Bridge | 91 | 180 / 180 | 0.116 (0.078) ms | **0.047 (0.016)** | 0.95 → 0.88 |
| Pyramid | 137 | 374 / 374 | 0.384 (0.143) | **0.285 (0.028)** | 1.75 → 1.64 |
| Soft Body | 193 | 3772 / 3772 | 0.790 (0.200) | **0.748 (0.108)** | 2.80 → 2.92（噪声内） |
| Pyramid XL | 1654 | 4868 / 4868 | 1.98 (1.18) | **1.40 (0.59)** | 8.70 → 8.35（噪声内） |

小场景上 LBVH 的宽相索引明显更快（哈希有固定层级/排序开销），大场景两者整步差距落在噪声里；默认保留哈希，LBVH 作为可切换的对照实现保留。

引擎的实际状态、测出来的数字和已知限制见 [`tools/baseline/README.md`](tools/baseline/README.md)（包括 529 刚体一步里各阶段的耗时、420 盒金字塔的穿透分布、以及一轮回归测试的原始输出）。

这不是一篇教程式的文章，能展示的就到这里。如果你也是做物理引擎的，你懂的。

[如果你真的很想知道的话](v2-6e8c45cac4eb0cacfc6b952955f36cae_1440w.jpg)

## 已知限制（Known limitations）

以下全部是**实测**结论，不是估计。数据、复现命令和原始输出在 [`tools/baseline/README.md`](tools/baseline/README.md)。

| 项目 | 现状 |
|---|---|
| 恢复系数 | **不模拟**。`body_set_param(BODY_PARAM_BOUNCE)` 只保存、只回读，碰撞完全非弹性 |
| 软体 | 自定义节点 `AVBDSoftBody3D` / `AVBDSoftWorld3D` 可用（刚体 + 关节格，见 `demo/scenes/soft_body.tscn`）；Godot 原生 `SoftBody3D` 的服务器侧 `_soft_body_*` 仍是 stub |
| Area3D | overlap 检测与 `body_entered` / `area_entered` 回调已实现（每步按形状做一次窄相，见 `_monitor_areas`）；RID / 参数 / layer / mask 可存取 |
| 每刚体形状 | **只支持一个** collision shape（多于一个 `WARN_PRINT` 后丢弃） |
| 形状类型 | box / sphere / cylinder 原生；capsule 用 cylinder 代理（半球帽变平）；凸包用它自己的 AABB 当盒子；concave / heightmap / world_boundary / separation_ray **排除出仿真**并告警 |
| 关节轴 | 只能沿**主轴**：非轴对齐的 hinge/slider/cone-twist 会被吸附到最近主轴（静默错误） |
| CCD | 没有。只有 speculative contact + `maxPenetration`=2 cm 兜底；5.3 m/s 落地实测最大瞬时穿透 **60 mm** |
| 静止穿透 | 420 盒双金字塔：p50 **8.0 mm**、p95 10.5 mm、p99 **11.3 mm**、最差稳定值 15.3 mm（盒子 1×1×0.5 m） |
| 冲击能量 | 5.3 m/s 落地的场景观测到 17.9 m/s 的弹飞速度：非弹性碰撞却在加速，冲击瞬间仍有过冲 |
| 刚性链预算 | 焊接/关节链每个 substep **至少 6 轮**才收敛；3 轮会自激发散（`tools/repro_joint.cpp`：140 m/s）。**总预算相同的情况下，substep 越多越不稳定** |
| 多线程 | 位级可复现（digest 在 1/2/3/4/8/12 线程下一致），但 **12 线程比 8 线程慢**（barrier 开销 > 并行收益），`physics/avbd/threads` 必须实测 |
| 宽相 | 默认分层哈希网格（`broad_phase_lbvh=false`），LBVH 可选；walk 内已做精确 AABB 过滤，候选对与 LBVH 相同。索引本身仍比 LBVH 慢（Pyramid XL 1.18 vs 0.59 ms），大场景整步差距落在噪声内，小场景约慢 8–10% |
| CI | 现在运行 `xmake run avbd_core_test`（138 checks）；README 顶部的 Test 徽章仍是编译徽章，不是测试徽章 |

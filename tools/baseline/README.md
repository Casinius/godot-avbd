# Baseline at `6ad6d63` (plan step 1)

Everything below was measured on this machine after making the tree build (see "build fix"
below). Nothing here is a fix: it is the reference point every later change must be
compared against.

```
host      : Linux x86_64, 12 cores
toolchain : g++ 16.2.0, xmake 3.1.2, -O3 -std=c++23, FP not fast-math
godot     : 4.7.2.stable.official
head      : 6ad6d63  ("fix collision: qmat tyz+-twx wrongly swapped")
params    : demo/project.godot as committed (240 Hz physics, 60 Hz internal tick,
            4 substeps x 3 effective rounds, beta_linear=1e5, projection on)
```

## Build fix (required, not cosmetic)

`6ad6d63` does not compile on gcc 16: `src/avbd/solver.cpp` uses `INT_MAX` without
including `<climits>` (older libstdc++ pulled it in transitively). Added the include.
The CI job only runs `xmake build` on a moving `ubuntu-latest`, so nothing in the repo
would have caught this.

**The CI never runs a single test.** `.github/workflows/c-cpp.yml` builds and runs
`xmake check`, which is not a target in `xmake.lua`. The README's "Test" badge is a
build badge.

## C++ core suite — `xmake run avbd_core_test`

```
132 checks, 0 failures          (tools/baseline/core_test.txt)
```

Digests recorded in that file (regenerate `tools/bridge_baseline.json` from them, its
current values are from an older build):

```
stack_300          = 0xfbac4f08cb1c977c
interleaved_300    = 0xbf1ce46c667164f3
mixed_240          = 0x521c838e09951796  (threads 1,2,3,4,8,12 all identical)
soft-body          = 0xbdf8ee1706e01184
```

Resting contact / stack penetration in the core suite (default params, 60 Hz):

```
rest contact penetration            0.00068 m
stack on ground (bottom box)        0.00108 m
stack interfaces (worst box-box)    0.00128 m
interleaved tight stack             0.00279 m
```

## Godot suites

| suite | result | file |
|---|---|---|
| `test_server.gd` | **30 checks, 0 failures** | `godot_test_server.txt` |
| `test_load.gd` | **PASS**, 626.8 ms per 400-box pyramid run (600 frames) | `godot_test_load.txt` |
| `test_solver.gd` | **38 checks, 1 failure** | `godot_test_solver.txt` |
| `test_constraints.gd` | **11 checks, 2 failures** | `godot_test_constraints.txt` |
| `test_pyramid.gd` | **9 checks, 3 failures** | `godot_test_pyramid.txt` |
| `test_penetration.gd` (30 s) | **19 checks, 2 failures** | `godot_test_penetration.txt` |

Failures, verbatim:

```
test_solver       hard joint error            max anchor separation = 8670.170898 m
test_constraints  over-constrained stays bounded   worst anchor error = 7024.7783 m
test_constraints  mixed chain joints hold          worst anchor separation = 908.7121 m
test_pyramid      no boxes below the ground        min y = -2281.433 m
test_pyramid      pile stays on the slab           max |xz| = 256.310 m
test_pyramid      centre of mass stays put         COM drift = 212.759 m in 10 s
test_penetration  no deep ground sink              min centre y = -2281.4326 m
test_penetration  penetration bounded              worst = 0.0665 m at t=9.5 s
```

What is actually measured, from `godot_test_penetration.txt`:

```
settled penetration  p50 0.0082  p95 0.0105  p99 0.0121  max 0.0132 m
                     (bodies still on the slab, 420-box double pyramid)
global max speed     206.1 m/s          bodies below the slab plane: 58 / 420
max |z|              255 m              max shift 2266 m
speed timeline       0.16 m/s @0 s -> 17 m/s @1.5 s -> 19 @10 s -> 203 @29.75 s
                     (a linear ramp: one body accelerated at ~8 m/s^2, never stopped)
```

So two separate things are true at once:

* the resting contacts are *better* than the doc comment in `test_penetration.gd` claims
  (that comment says p50 1.4 cm / p99 2.0 cm; measured now: 0.82 / 1.21 cm);
* a handful of boxes are **ejected** during the initial drop (17 m/s by t=1.5 s, against
  a 5.3 m/s impact), leave the 80x20 m slab, and then free-fall for the rest of the run.
  Every "below the ground" / "on the slab" / "COM drift" failure is that, not a contact
  problem.

## The same suites before `6ad6d63`

`6ad6d63` fixed `qmat()` (a matrix that was not a rotation) and `orthonormal()` (a frame
whose rows were its columns). Both fixes are correct. Their *effect* on the suites is not
uniformly positive, so the parameter set was clearly tuned against the broken build:

| suite | `8d1e4fc` (before) | `6ad6d63` (now) |
|---|---|---|
| `avbd_core_test` | 9 failures (every sphere/cylinder test flew off) | **0 failures** |
| `test_solver` | 2 failures: stack height 0.09485 m, stack drift 0.03305 m | 1 failure: hard joint 8670 m |
| `test_constraints` | **0 failures** | 2 failures: joint chains explode |
| `test_pyramid` | **0 failures** (max speed 0.006 m/s, zero drift) | 3 failures (206 m/s ejections) |

Built with `git worktree add /tmp/avbd-pre 8d1e4fc` + the same `<climits>` fix. The
pre-fix build has no `get_solver_max_penetration` diagnostic (added in `6ad6d63`), so
`test_penetration` cannot run on it.

## Root cause of the ejections: the dual variable has no bound

Isolated in pure C++ (`tools/repro_joint.cpp`, no Godot): a 4-link chain of fully locked
`GenericJoint` welds under a static anchor, at the demo's substep schedule.

```
g++ -O2 -std=c++23 -Isrc -I<thread-pool>/include tools/repro_joint.cpp src/avbd/*.cpp \
    -o /tmp/repro_joint -pthread
/tmp/repro_joint <substeps> <rounds> <links> <lockAngular> <dtDiv> <betaLin>
```

| substeps | rounds/substep | dt | result |
|---|---|---|---|
| 1 | 12 | 1/240 | stable (vz ~0.01) |
| 4 | 12 | 1/240 | stable |
| 4 | **3** | 1/240 | **explodes**: vz -> 1400 m/s, z -> -4500 m |
| 4 | 3 | 1/960 | stable |
| 4 | 3 | 1/120 | explodes |

and, at 3 rounds, by `betaLin`:

```
betaLin=100     stable (k saturates at 675, chain sags 4 cm)
betaLin=1000    explodes
betaLin=10000   explodes
betaLin=100000  explodes        <- the demo's setting
```

The trace shows the mechanism (`lambda` of the anchor joint's Z axis):

```
t=3    kZ=  3600  lamZ=     80   C= 0.0046
t=20   kZ=  3600  lamZ=  -1416   C= 0.0434
t=60   kZ=  3600  lamZ=-1.96e7   C=-221.9
```

`penalty` is capped by `GenericJoint::penaltyLimit` (3600 N/m here) but `lambda` is not
capped by anything: `lambda = penalty*C + lambda` accumulates the residual `C` every
round. When the primal cannot drive `C` to zero inside the round budget, `lambda` grows
without bound, the next step's force is enormous and signed against the motion, and the
chain is thrown. `ldlt6_solve`'s "unsane update" guard (reject if > 10 m) is three orders
of magnitude too coarse to catch it.

Same failure mode, same signature, in the pyramid: with
`--sweep-config=rounds-8` (32 rounds/tick instead of 12) the pyramid goes from

```
global_max_speed 206.1 -> 18.7 m/s      bodies_below_slab 58 -> 0
p50 / p99        0.0082/0.0121 -> 0.0070/0.0080 m
```

i.e. more iterations mask it, which is exactly what an unbounded dual does. Contacts are
protected by `maxPenetration = 0.02` clamping the *error*; joints have no equivalent and
contacts have no bound on `lambda` either, only on `C`.

## Consequence for the plan

Step 2 of the plan (`contactOnlyOverflow = false`) is not the top priority any more.
The order should be:

1. bound the dual (`lambda`) and/or the per-step joint/contact force, and make the
   `ldlt6` sanity bound scale with `dt` and the body's mass instead of a fixed 10 m;
2. re-tune `betaLin` / rounds against the *fixed* build;
3. only then A/B `contactOnlyOverflow`, the projection pass, and `(1-alpha)`.

## Reproducing this file

```sh
xmake f -m release -y && xmake build avbd && xmake build avbd_core_test
xmake run avbd_core_test > tools/baseline/core_test.txt
for s in test_server test_load test_solver test_constraints test_pyramid; do
    godot --headless --path demo --script res://tests/$s.gd > tools/baseline/godot_$s.txt
done
godot --headless --path demo --script res://tests/test_penetration.gd -- --pen-seconds=30 \
    > tools/baseline/godot_test_penetration.txt
```

`demo/tests/probe_chain.gd` prints the per-tick state of the welded chain through the
Godot path (it is what showed that the chain rings, then accelerates, and that the
bodies are *not* asleep while it happens).

---

# Step 2/3 result (this working tree)

## The stability knob is rounds-per-substep, not the penalty cap

At a **fixed 12 rounds per tick**, whether a 4-link welded chain settles depends only on how
those rounds are split across substeps (`tools/repro_joint.cpp`, betaLin 1e5, dt 1/240):

```
substeps x rounds   |vz|max after 300 ticks
1 x 12              0.07 m/s     stable
2 x 6               0.06 m/s     stable
3 x 4               0.31 m/s     marginal
4 x 3             140.49 m/s     diverges      <- the demo's schedule before this change
6 x 2             413.58 m/s     diverges
12 x 1           1310.33 m/s     diverges
```

Fewer substeps = more rounds inside one linearisation, and the augmented-Lagrangian dual
converges *inside* a linearisation. Interleaved substepping re-linearises the constraint
(and rebases it) every substep while the dual keeps accumulating, so below ~6 rounds per
substep the dual never converges and the constraint becomes a launcher.

Capping `GenericJoint::penaltyLimit` was tried first (swept over 3/6/8/12 rounds and caps
0.25/0.125/0.0625/0.032 of `omega*dt`): the stable region is not monotone in the cap, and
every cap low enough to be stable at 3 rounds sags a weld by >3 cm. That is a tuning
dead end; the round budget is the actual fix.

**Change:** `demo/project.godot` runs `substeps=2` instead of 4 (the automatic split then
gives 6 rounds/substep, i.e. the same 12 rounds/tick, same cost), and the no-op
`rounds_per_substep=2` is gone. `Solver::stepInterleaved` now treats `roundsPerSubstep` as
an override rather than a floor, so the setting can no longer be silently overruled.

Result at equal cost (626 ms per 400-box/600-frame pyramid run, unchanged):

| | before (4x3) | after (2x6) |
|---|---|---|
| `test_pyramid` | 3 failures, 206 m/s ejections, 58 bodies off the slab | **0 failures** |
| `test_solver` | 1 failure (joint 8670 m) | **0 failures** |
| `test_constraints` | 2 failures (chains 7024 / 908 m) | **0 failures** |
| penetration p50 / p99 | 8.2 / 12.1 mm | 8.0 / 11.3 mm |
| worst *settled* penetration | 15.3 mm | 15.3 mm |
| `test_penetration` | 2 failures | 1 failure: 60 mm transient at the 5.3 m/s drop impact |

## `contactOnlyOverflow` A/B (kept at `true`)

```
                    p50 pen   p99 pen   worst settled   max speed   cost
true  (Jacobi)       8.0 mm    11.3 mm      15.3 mm       17.9 m/s   626 ms
false (colour GS)    7.4 mm     8.8 mm      36.3 mm       18.4 m/s   626 ms
```

Gauss-Seidel has the better median but a systematic tail (a pair the colour order visits
last settles 3.6 cm deep: the update order is a property of the scene, not of the load).
Jacobi averages that away. Kept `true`; exposed as `physics/avbd/contact_gs`.

## Bugs fixed in this step (with regressions)

* `Solver::pick` tested every body as a box of `size*0.5`: a sphere is stored as
  `{radius,0,0}`, so it became a zero-thickness slab (rays could not hit one) and a
  cylinder's radius was halved. Now per shape, with an exact ray/sphere.
  Tests: `pick hits a sphere off its centre`, `pick lands on the sphere surface`,
  `pick hits a cylinder off its axis`, `pick misses outside the cylinder radius`.
* The broad phase emitted candidates for static-static pairs; `ManifoldStore::frozen()`
  then kept those zero-contact manifolds alive forever (one slot + 8 contact records per
  pair). Test: `static-static pairs make no manifold` (was 5, now 0).
* The AABB wake-up wrote `Rigid::sleeping`/`stillFrames` from the *parallel* pair filter:
  a data race, and the order-dependent part was invisible because it was "idempotent".
  Moved into the serial manifold-creation loop.
* `Joint::updateDual` ramped `penaltyAng` twice per pass, once inside the *linear* block
  using the linear error (metres) as an angular one and `betaAng`. Inherited from the
  reference; the angular block is now the only writer.

All four are caught by the new checks: reverting `src/avbd/solver.cpp` makes 4 of them fail.

```
avbd_core_test   138 checks, 0 failures   (132 before, +6 new)
test_pyramid     9/9      test_solver  38/38     test_constraints 11/11
test_server     30/30     test_load    PASS
test_penetration 19 checks, 1 failure (60 mm transient at impact only)
```

---

# Step 4 result: where a step actually goes

`Solver::lastStepTimes` now records the four phases (std::chrono, 4 reads per step) and the
`get_broadPhase_time()` / `get_colourGraph_time()` / `get_solve_time()` /
`get_finish_time()` accessors return them instead of a hardcoded 0. The core bench prints
the breakdown of its last step:

```sh
xmake run avbd_core_test --threads N --bench
```

529 boxes, 5471 contacts, 8 iterations, 60 Hz, per-step milliseconds:

| threads | broad | graph | solve | finish | total |
|---|---|---|---|---|---|
| 1 | 1.27 | 0.85 | 5.51 | 0.81 | 8.07 |
| 2 | 1.15 | 0.39 | 3.39 | 0.81 | 5.50 |
| 4 | 1.22 | 0.37 | 2.30 | 0.81 | 4.37 |
| 8 | 1.38 | 0.32 | 1.64 | 0.87 | 4.17 |
| 12 | 1.75 | 0.32 | 1.54 | 1.05 | 4.59 |

Three things this says:

* The iteration region is the cost (5.5 -> 1.5 ms, 3.6x on 12 threads - it is barrier
  bound, not bandwidth bound) but the **broad phase does not scale at all** (1.27 -> 1.75
  ms) and is 16-38% of a step. That is where an LBVH-vs-hash decision actually matters, and
  it is the biggest serial block left.
* `finish` costs ~0.8-1.0 ms, of which most is `projectContacts` - the second (Jacobi
  position) solver that runs every step. That is the price of step "position projection on
  by default" (see the review), and it is comparable to the whole broad phase.
* 12 threads is *slower* than 8 on this scene: the extra barriers cost more than the added
  parallelism. `physics/avbd/threads` should be benchmarked, not maxed.

## The hash-vs-LBVH A/B, finally isolated

The old claim ("the spatial hash is faster than the dynamic LBVH") compared this branch
against `main`, which is a different engine. So the LBVH was ported: `src/avbd/bvh_broadphase.*`,
from `main:src/avbd/bvh/`, changed only in what a comparison needs - Eigen -> linalg, and
**leaves sized from the same swept world AABBs the hash receives** (`half`/`centre` spans)
instead of `body->radius` bounding spheres. Both broad phases now emit the same packed pair
candidates into the same `sweepPairs`, and everything downstream (radix sort, unique, exact
AABB/layer/constraint filter, manifold slot creation) is shared.

`Solver::broadPhaseMode` selects one; `--bvh` in the core test does the same from the CLI.

**The simulation is bit-identical either way** - `stack_300 = 0xfbac4f08cb1c977c`,
`interleaved_300 = 0xb03adeb37742850a`, `mixed_240 = 0xef9ebb5b27a25c27` with both - and the
contact count matches on every scene, so the two find exactly the same real pairs.

135+ bodies, 8 threads, per-step milliseconds (`--bench-scene NAME`):

| scene | bodies | pairs hash / lbvh | broad hash (index) | broad lbvh (index) | step hash -> lbvh |
|---|---|---|---|---|---|
| Stack | 11 | 25 / 10 | 0.031 (0.020) | **0.010 (0.000)** | 0.108 -> 0.091 |
| Bridge | 91 | 1295 / 180 | 0.307 (0.090) | **0.049 (0.017)** | 1.145 -> 0.841 |
| Pyramid | 137 | 2865 / 374 | 0.515 (0.200) | **0.201 (0.029)** | 1.611 -> 1.345 |
| Soft Body | 193 | 12099 / 3772 | 1.067 (0.169) | **0.623 (0.098)** | 2.518 -> 2.242 |
| Pyramid XL | 1654 | 45603 / 4868 | 3.366 (1.185) | **0.953 (0.582)** | 10.019 -> 8.243 |

Contacts after 140 steps are identical in every row (1438 / 19447 / 40 / 293 / 3800). The
LBVH wins **despite being serial** while the hash's per-body work and walk run on 8 threads.

The reason is a defect in the hash, not a strength of the tree: it walks the 27-cell
neighbourhood at *every* level from the body's own up to the coarsest. A scene with one huge
body - the ground - puts every other body's coarse-level walk in a ground-sized cell, so the
candidate list inflates 6-9x (45603 vs 4868 on Pyramid XL) and the exact filter plus the pair
sort pay for all of it. Fixing the walk (own level only, with coarse bodies resolved through
the fine cell's parent) is the way back to the grid.

The default is now the LBVH (`broadPhaseMode = 1`, `physics/avbd/broad_phase_lbvh`). At the
demo's own scale (420 boxes) the difference is inside the noise: `test_load` 629.9 ms vs 626.8,
`probe_engines` 1.43 ms/frame either way, penetration p50/p99 unchanged at 8.08/11.30 mm.

> Later cleanup: the hashed grid was kept and optimized instead of removed. The walk now
> emits a same-level pair once, skips empty levels, and tests the swept AABBs before
> writing the candidate, so the candidate list matches LBVH exactly (4868 vs 45603 on
> Pyramid XL). `broadPhaseMode` remains, with the hash as the default (`0`) and LBVH as
> the `1` alternative. See README.md for the post-optimization numbers.

# Step 5 result: the `(1-alpha)` contact damping is load-bearing

`ManifoldStore::constraintError` returns `C = (1-alpha) * C0 + B * dq`, so at the default
alpha = 0.99 only 1% of the step-start overlap is visible to the constraint. Swept on the
420-box pyramid (15 s, `physics/avbd/alpha`), `p50`/`p99` = per-contact geometric overlap:

| alpha | p50 | p99 | worst settled | max speed | note |
|---|---|---|---|---|---|
| 0.0 | 4.0 mm | 5.4 mm | 5.5 mm | **96 m/s** | full AL, best penetration, launches boxes |
| 0.9 | 5.2 mm | 6.8 mm | 7.4 mm | **111 m/s** | |
| **0.99** | **8.0 mm** | **10.3 mm** | **10.6 mm** | **17.9 m/s** | default |
| 1.0 | 19.3 mm | 22.8 mm | 23.8 mm | 79 m/s | pure bracket: no positional term at all |

So it is not a softening bug that can simply be removed: `alpha = 1` (removing it) makes both
penetration *and* stability worse, and lowering it buys penetration at 5-6x the impact
launch speed. The default sits on a stability optimum. It is inherited from the reference
implementation (`main:src/avbd/manifold.cpp:212` has the same expression), and the same
`alpha` also drives the dual warm-start decay - one parameter in two roles, which is why it
cannot be tuned in one direction.

# Step 6 result: the block inertia is now world-frame

`Block` used `diag(moment)/dt^2` (the body-frame inertia) against world-frame constraint
Jacobians and a world-frame angular update. It now builds `I_w = R diag(moment) R^T` and
uses `I_w/dt^2` for both the Hessian and the right-hand side, matching what
`projectContacts`' `angularCompliance`/`angularResponse` already did.

Measured (no regression, so it stays):

```
core suite            138 checks, 0 failures      digests changed (expected)
bench (529 bodies)    solve 1.646 ms/step  (1.642 before)   contacts 5416 (5471 before)
test_pyramid 9/9   test_solver 38/38   test_constraints 11/11
penetration (15 s)    p50 8.08 mm  p99 11.30 mm  max speed 18.4 m/s
test_load             627 ms per 400-box pyramid run (626.8 before)
```

# Step 7 result: dead code and stale docs

Verified first (the earlier review guessed wrong here): `newtonRatio`, `stiffnessDecay` and
`adaptiveBeta` are **live** - the round split, the dual-tail stiffness decay and the
per-contact beta multiplier all still run, they are just off/default by default. Nothing
deleted there.

Actually removed (all unreferenced, `grep` count 1 = the definition only):

* `Solver::computeCollisionMargin`, `Solver::currentCollisionMargin`,
  `COLLISION_MARGIN_DAMPING`, `COLLISION_MARGIN_MAX` - a dynamic-margin idea that was never
  wired to anything; `currentCollisionMargin` was assigned nowhere.
* `ShapeType::Cylinder`'s doc block said the axis is local +Y and the height is `size.y`;
  every implementation (Rigid, aabbHalfExtent, makeShape, the server) uses **+Z / size.z**.
* `ManifoldStore` honoured `STICK_DRIFT_HEAL` directly instead of `Solver::stickDriftHeal`,
  so the documented per-solver knob did nothing. Wired up (default unchanged).
* `Shape::axis` had no default, so five `Shape{...}` aggregates in `test/core_test.cpp`
  warned under `-Wextra` - the tree was not warning-clean despite `xmake.lua` claiming it is.
  Now defaulted; the core test builds silently.
* README: the "dynamic BVH / SAH refit" claim was false (the broad phase is a hashed grid),
  and "no need to worry about maturity" is replaced by the measured limitations table.

# Step 8 result: the README now states what is not supported

`## 已知限制（Known limitations）` in README.md: no restitution, soft bodies and Area
monitoring stubbed, one shape per body, capsule/hull proxies, joint axes snapping to the
nearest principal axis, no CCD (with the measured 60 mm impact transient), the measured
penetration distribution, the rigid-chain round budget, the 12-vs-8-thread anomaly, the
non-scaling broad phase, and that CI never runs a test.

# Step 9 result: same scene, same frames, AVBD vs Jolt

`demo/tests/probe_engines.gd` loads the same `pyramid_dbl.tscn` (420 boxes), disables
sleeping, runs the same frame count, and reports wall time plus the lowest box's sink below
its authored resting height. Both engines run at the project's 240 Hz physics rate; AVBD
additionally runs 60 Hz internal ticks with 2 substeps x 6 rounds (12 rounds/tick, 720
rounds/s). Jolt runs its own default iteration counts.

| frames (sim time) | engine | ms/frame | ground sink | max speed |
|---|---|---|---|---|
| 300 (5 s) | AVBD | 1.48 | 10.33 mm | 14.2 m/s |
| 300 (5 s) | Jolt | 1.23 | 0.87 mm | 12.8 m/s |
| 1200 (20 s) | AVBD | 1.42 | **2.95 mm** | 10.5 m/s |
| 1200 (20 s) | Jolt | 1.31 | **3.13 mm** | 10.1 m/s |

At 20 s the two are within 8% on wall time and AVBD is 0.2 mm *better* on ground sink. The
5 s row is the interesting one: Jolt holds 0.9 mm through the initial impact while AVBD is at
10 mm, i.e. AVBD's steady state is fine and its weakness is the transient (same story as the
60 mm impact spike in step 1). Both piles are still moving at 20 s (the double pyramid's
upper row is unsupported), so neither number is a settled-pile measurement.

`tools/compare_engines.py` (wall time over 300 fixed frames, minus an empty-scene control)
agrees on the order of magnitude: `pyramid_dbl` AVBD 8.50 ms vs Jolt 7.79 ms per frame,
`pyramid` 6.17 vs 6.56, `stack` 2.32 vs 2.63 - absolute numbers there include the demo's
per-frame driver work, which is why they are 5x the probe's.

## Step 9: the three-engine comparison, same scene and frames

`demo/tests/probe_engines.gd` loads `pyramid_dbl.tscn` (420 boxes), disables sleeping, runs
1200 physics frames (20 s of sim) and reports wall time plus the lowest box's sink below its
authored resting height. Only `physics/3d/physics_engine` changes between runs.

| engine | ms/frame | ground sink | max speed | off_slab |
|---|---|---|---|---|
| GodotPhysics3D | 2.19 | 11.82 mm | 10.19 m/s | 0 |
| Jolt Physics | **1.38** | 3.13 mm | 10.08 m/s | 0 |
| **AVBD** | 1.45 | **2.95 mm** | 10.53 m/s | 0 |

AVBD is 5% behind Jolt on wall time, 51% ahead of GodotPhysics3D, and has the lowest ground
sink of the three. Budgets are not equal: AVBD runs 240 Hz physics with 60 Hz internal ticks
and 2 substeps x 6 rounds (12 rounds/tick, 720 rounds/s of 6x6 block solves); Jolt runs its
own default iteration counts at 240 Hz; GodotPhysics3D runs sequential impulses at 240 Hz
with `3d/solver/solver_iterations=4`.

Caveats, stated because the earlier 5 s sample of this same run says the opposite: at 300
frames (5 s) AVBD showed 10.3 mm sink against Jolt's 0.87 mm, i.e. AVBD's steady state is
fine and its weakness is the impact transient (the same story as the 60 mm spike in
`test_penetration`). Both piles are still moving at 20 s - the double pyramid's upper row is
unsupported - so none of these numbers is a settled-pile measurement.

`tools/compare_engines.py` (wall time over 300 fixed frames minus an empty-scene control)
agrees on the order of magnitude: `pyramid_dbl` AVBD 8.50 vs Jolt 7.79 ms per frame,
`pyramid` 6.17 vs 6.56, `stack` 2.32 vs 2.63. Its absolute numbers include the demo's
per-frame driver work, which is why they are several times the probe's.

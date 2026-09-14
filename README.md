# godot-avbd

AVBD (Augmented Vertex Block Descent) physics for Godot 4, shipped as a GDExtension
and built with **xmake** against the `godotcpp4` package from xmake-repo.

AVBD is the augmented-Lagrangian formulation of Vertex Block Descent (Giles, Diaz,
Yuksel, SIGGRAPH 2025): an unconditionally stable, parallelisable position-based
solver that handles hard constraints at infinite stiffness, stiff systems with large
stiffness ratios, and friction contacts. The solver core here is a direct port of the
authors' reference implementation (`savant117/avbd-demo3d`, MIT); the Godot node
layer, build system, tests and demo are new.

```
src/avbd/     solver core, no Godot dependency: rigid bodies, contacts
              (box-box SAT + face clipping), joints, springs, collision filtering,
              the constraint-graph colouring that drives the parallel update, and
              job_pool.hpp (the adapter over BS::thread_pool)
src/nodes/    Godot node classes + the Y-up <-> Z-up conversions
tools/        the reference's test scenes, usable without Godot
test/         numeric assertions for the core
demo/         Godot 4 project: five scenes, a demo driver, headless tests
```

## Requirements

| | |
|---|---|
| xmake | 3.x, with the `godotcpp4` (builds godot-cpp via `scons` on first use; 4.1 is what this project pins) and `thread-pool` packages available from xmake-repo |
| Godot | 4.1 or newer. Verified against 4.4.1 and 4.6.3 |
| compiler | any C++20 compiler (verified with GCC 16) |

The solver's parallelism uses [BS::thread_pool](https://github.com/bshoshany/thread-pool)
(MIT, header-only; `add_requires("thread-pool v5.1.0")`). It provides the job queue, the
worker lifetime and the fork-join handshake; the solver supplies only the part a library
cannot know, which is which bodies may be updated at the same time.

## Build

```sh
xmake f -m release && xmake     # libavbd.so, copied to demo/bin/
xmake run avbd_core_test        # solver-only tests: no engine, no window
```

Both targets build clean at `-Wall -Wextra -Wshadow -Wnon-virtual-dtor`. Sanitizer builds
are configuration flags, so a clean tree can be checked with one command each:

```sh
xmake f -m debug --asan=y --ubsan=y && xmake    # address + undefined behaviour
xmake f -m debug --tsan=y && xmake              # data races (needs its own build dir)
```

`xmake` builds two targets:

* `avbd` - the GDExtension shared library, linked against the static `godotcpp4`
  build. The post-build step copies it to `demo/bin/libavbd.linux.x86_64.so`, which
  is what `demo/avbd.gdextension` points at.
* `avbd_core_test` - a standalone binary: numeric solver tests, scene runner, and a
  528-box benchmark.

## Run the demo

```sh
Godot --headless --import --path demo        # once per checkout: registers the extension
Godot --path demo                            # stack scene
Godot --path demo res://scenes/pyramid.tscn  # also: friction, soft_body, bridge
```

Controls: **left click** punches along the camera ray, **Space** pauses, **R** resets
every body to its scene pose. For unattended runs:

```sh
Godot --headless --path demo res://scenes/soft_body.tscn -- --steps=300
```

## Tests

```sh
xmake run avbd_core_test                                       # 41 checks, solver only
Godot --headless --path demo --script res://tests/test_load.gd # classes + property round-trip
Godot --headless --path demo --script res://tests/test_solver.gd       # 41 checks, one constraint at a time
Godot --headless --path demo --script res://tests/test_constraints.gd  # 74 checks, constraints composed
```

The core binary also has `--bench`, `--bench-sweep` (ms/step per thread count), `--colour`
(colour count and widest group per scene), and `--threads N`, `--scene NAME`, `--steps N`.

`test_solver.gd` covers single features (rest contact, stacking, friction, joints,
raycasts, impulses, soft bodies). `test_constraints.gd` covers *combinations* - rigid
joints plus springs plus springs on the same bodies, chains that mix hard and compliant
links, over-constrained pins, collision suppression, fracture under load, a plank bridge,
a lattice on a spring mattress, a teleported static driver, `substeps`, and that the thread
count does not change a result.

Both suites share `demo/tests/avbd_harness.gd` and accept:

| flag | meaning |
|---|---|
| `--list` | print scenario names |
| `--only=NAME` | run one scenario in isolation |
| `--digest` | run each scenario twice and compare its evidence, plus a whole-flow digest |
| `--strict` | stop at the first failure |

Each scenario prints the measurements it made (`evidence: {...}`), so a run is a record
of the numbers, not just pass/fail. Where a result is predictable from theory the
assertion uses the prediction rather than a fitted bound - a spring mattress must sag by
`m*g/(n*k)`, a joint-plus-spring pair must hold its link length to millimetres, a spring
chain link must stretch by `load/k`, and an impulse of `J` on a mass `m` must give `J/m`.

`demo/tests/capture.gd` renders a scene to a PNG for visual checks (needs a real
rendering device, so do not pass `--headless`):

```sh
Godot --path demo --script res://tests/capture.gd -- \
    --scene=res://scenes/stack.tscn --frames=300 --out=/tmp/stack.png
```

## API

### AVBDWorld3D

Root of a simulation: owns the solver and steps it once per physics tick. Every
`AVBDRigidBody3D`, `AVBDJoint3D`, `AVBDSpring3D`, `AVBDIgnoreCollision3D` and
`AVBDSoftBody3D` in its subtree is picked up automatically, in scene order. This
node's transform defines the simulation's coordinate system; the rest of the scene
is unaffected by it.

| property | default | meaning |
|---|---|---|
| `gravity` | 9.8 | m/s² along -Y |
| `iterations` | 10 | solver iterations per step |
| `alpha` | 0.99 | constraint stabilisation; higher = smoother, less energetic |
| `beta_linear` / `beta_angular` | 10000 / 100 | penalty ramping, i.e. how fast contacts and joints stiffen |
| `gamma` | 0.999 | warm-start decay of the dual variables |
| `substeps` | 1 | solver steps per physics tick (each of `delta / substeps`) |
| `threads` | 0 | worker threads: 0 = one per hardware thread, 1 = everything on the calling thread, N = N workers. Does not change the result |
| `paused` | false | stop stepping (bodies stay where they are) |

| method | |
|---|---|
| `rebuild()` | re-create the solver from the scene graph, keeping current body poses. Called automatically when the subtree changes |
| `raycast(origin, direction, max_distance = 10000)` | closest hit, or an empty Dictionary. Keys: `body`, `position`, `distance` (plus `cell` for lattice soft bodies). Static bodies are hit too |
| `get_body_count()`, `get_force_count()`, `get_contact_count()`, `get_contact_point_count()`, `get_step_time_usec()`, `get_step_count()`, `get_thread_count()` | statistics |

### AVBDRigidBody3D

A box-shaped body. Its transform when it joins the simulation is its initial pose;
afterwards the solver owns a dynamic body and writes the result back every tick.

| property | default | meaning |
|---|---|---|
| `size` | (1, 1, 1) | box extents in metres (full widths) |
| `density` | 1.0 | kg/m³; mass and inertia are derived from `size` |
| `friction` | 0.5 | two bodies in contact use the geometric mean |
| `static_body` | false | immovable collider. Its node transform stays authoritative, so moving one in the editor or by script moves the collider |
| `initial_velocity` / `initial_angular_velocity` | 0 | applied when the body enters the simulation |
| `linear_velocity` / `angular_velocity` | - | current velocity, world axes |
| `visualize_shape`, `shape_color` | false | debug box mesh, created at runtime only |

Methods: `get_mass()`, `apply_impulse(impulse, position_offset = 0)`,
`teleport(position, rotation = identity)`. `teleport` moves a dynamic body and stops
it; opening a scene away from the sim while the bodies are already moving, e.g. for a
respawn.

### AVBDJoint3D

Ball-and-socket constraint, optionally angular and breakable.

| property | default | meaning |
|---|---|---|
| `node_a`, `node_b` | empty | the constrained bodies, as NodePaths relative to the joint (or absolute) |
| `anchor_a`, `anchor_b` | 0 | attachment points. `anchor_a` is a body-local offset, **or a world-space position when `node_a` is empty**, which pins body B to a point in space |
| `linear_stiffness` | -1 | N/m. **negative = infinitely stiff (hard)**, 0 disables the linear constraint, positive gives a soft joint with its penalty capped there |
| `angular_stiffness` | 0 | negative = also locks rotation, 0 = free rotation, positive = soft angular spring |
| `fracture_force` | -1 | angular force above which the joint breaks; negative = unbreakable |
| `broken` | false | read-only in practice: true once the joint has broken. Breaking is permanent for the node |

A joint that breaks is deleted by the solver, which also means its two bodies start
colliding with each other again. Because the solver forgets deleted joints, the node
remembers: `is_broken()` stays true across a `rebuild()` (adding or removing any body),
and a broken joint is recreated without a solver force instead of silently re-attaching.
`break_joint()` breaks it on demand and `broken = false` mends it for the next rebuild.

Every constraint node also exposes `is_simulated()`, which reports whether the world
managed to create its solver force - if a `NodePath` is wrong, or the endpoint body is
outside the world, the constraint is skipped and this is false.

### AVBDSpring3D

Linear spring between two bodies (both endpoints required - anchor one to a static
body to hang a spring from the world). `anchor_a` / `anchor_b` are body-local
offsets, `stiffness` is in N/m, and `rest_length` of -1 derives the rest length from
the bodies' pose when the spring is created.

Anchors that coincide at spawn with `rest_length = 0` give a "face spring" that holds
two touching faces together with a finite stiffness: under a load `F` the faces separate
by exactly `F / k`, which is what the `mixed_chain` scenario checks.

### AVBDIgnoreCollision3D

Suppresses contacts between two bodies. Bodies joined by *any* constraint already
skip each other, so this is for pairs that would otherwise collide spuriously.

### AVBDSoftBody3D

A soft body in the formulation the AVBD authors use: a `dimensions` grid of jointed
boxes, with diagonal neighbours' collisions suppressed. The cells are real solver
bodies, so the lattice collides with the rest of the scene. It draws through a single
`MultiMesh`, so a 4x4x4 blob is one node and one draw call. The lattice is centred on
the node's origin.

| property | default | meaning |
|---|---|---|
| `dimensions` | (4, 4, 4) | lattice resolution |
| `box_size` | (0.25, 0.25, 0.25) | cell extents |
| `gap` | 0 | extra space between cells (0 = touching) |
| `density`, `friction` | 1.0, 0.5 | per cell |
| `stiffness_linear`, `stiffness_angular` | 1000, 250 | the links between cells |

`get_cell_count()` and `get_cell_transform(index)` expose individual cells in the
node's local space.

## Conventions

* **Coordinates.** Godot is Y-up; the solver core is Z-up with gravity along -Z. The
  mapping is a +90° rotation about X (`sim = (x, -z, y)`), applied in exactly one
  place, `src/nodes/godot_convert.hpp`. Note that this flips one sign: Godot +Z runs
  along solver -Y, so any offset expressed in solver space has to go through the same
  conversion (this was a real bug in the soft-body lattice anchors).
* **Units.** Metres, kilograms, seconds, radians. The solver parameters (`beta_*`,
  `alpha`, `gamma`) are scale-dependent, exactly as the reference implementation warns.
* **Determinism.** Single-threaded *update order* with no randomness: the same scene stepped
  the same number of times gives bit-identical results - within a process, across processes
  (checked for the core binary and through the engine), and for every `threads` setting.
  The suites print digests to compare between runs; the current 300-step stack digest is
  `0xabe38013f0c29277`.
* **Parallelism.** Bodies are grouped ("coloured") so that no two bodies in a group share a
  constraint; a group is then updated concurrently, one group at a time. That is exactly the
  condition for the updates to be independent, which is why increasing `threads` cannot
  change a single number - only how long it takes. The colouring fixes the update order, so
  a scene gives the same trajectory on one thread and on twelve. How much of a scene can run
  at once is decided by the widest group, not the thread count: a settled 10-box stack has 2
  groups (a stack is a path), a 528-box pyramid has 4 whose widest holds half the boxes, so
  ~2x is its ceiling.

  | phase | |
  |---|---|
  | inertial/warm-start, primal update, velocity update (BDF1) | parallel, by body |
  | force warm-start and dual update | parallel, by force (deletion stays in one thread) |
  | broad-phase pair test and manifold creation | serial - an O(n²) scan and order-sensitive list splicing |
  | colouring | serial, ~0.01 ms for 500 bodies |

  Measured on a 12-core machine, 528-box pyramid, 10 iterations: 9.0 ms/step on one thread
  and 3.95 ms on eight or twelve (2.3x).
* **Subtree changes.** The world scans its subtree each tick and rebuilds the solver
  when the set of nodes changed, carrying over the current pose of bodies that stay.
  Adding or removing a body at runtime therefore does not reset the simulation.
* **Constraint counting.** `get_force_count()` counts everything the solver holds,
  including contact manifolds. To count only declared constraints, subtract
  `get_contact_count()`; that is what the test suites' pre-flight does.
* **Stabilised constraints do not snap.** Teleporting a static body that a joint is
  attached to leaves an instantaneous violation which the solver closes over seconds
  (`alpha` controls the rate), rather than moving the other body in one step. The
  `static_driver` scenario measures this: a 4 m teleport is fully recovered in 437
  steps (7.3 s) with the hitch accurate to 0.6 mm afterwards.
* **Free joints are pendulums.** With `angular_stiffness = 0` a joint only pins a point,
  so the body hangs and swings about it. The `static_driver` scenario sets
  `angular_stiffness = -1` to get a rigidly hitched body that tracks its driver exactly.

## Measured behaviour

Numbers from this machine (i7-9750H, 12 threads, release build), `dt = 1/60`:

| scene | bodies | contacts | ms/step |
|---|---|---|---|
| stack (11 boxes) | 11 | 40 points | 0.12 |
| friction (ramp + 8 boxes) | 10 | 33 points | 0.07 |
| bridge (20 planks + loads) | 29 | 46 points | 0.27 |
| pyramid (16 rows) | 137 | ~1030 points | 2.1 |
| soft-body demo (3 lattices, 14 iterations) | 193 | ~340 points | 5.5 |
| benchmark, 528 boxes, 10 iterations | 528 | 4545 points | 3.95 (9.0 on one thread) |

The same 528-box benchmark before the parallel update and the `constrainedTo` fix measured
53 ms/step; it is now 3.95 ms, of which the parallel part is 2.3x faster than its serial
form.

Accuracy, as asserted by the tests:

* a box dropped on a slab rests with ~9 mm penetration (the model is a soft penalty
  with a 1 cm collision margin); stacking compresses ~2 cm per interface, and the
  stack does not drift sideways at all.
* stiction: on a 20° ramp with mu = 0.5 the box creeps ~0.7 mm per half second, versus
  4.7 m per half second at mu = 0.05.
* hard joints (infinite stiffness, `anchor` on both sides): ~2 mm separation error on a
  4-link chain, ~1 mm for a body pinned to a world point.
* an impulse of 10 N·s on a 1 kg box yields exactly 10 m/s, and at its top face an
  angular velocity of -30 rad/s, i.e. the inertia tensor is used correctly.
* contacts are penalty-based, so an *initial* constraint violation decays over a second
  or two rather than snapping; spawn bodies at their constrained pose.

Constraint combinations (`test_constraints.gd`):

| composition | measured |
|---|---|
| rigid joint + 400 N/m spring on one body | link length held to 0.18 mm while the spring pulls 277 N; bob swings 67°; settles to 0.15 mm/s |
| rigid link + spring link chain (k = 200 N/m) | rigid spans 0.99998 m (1 m nominal); spring spans 1.147 / 1.049 m for loads of 3 / 1 kg (`load/k` = 1.147 / 1.049) |
| one body on two world joints | consistent anchors: 1 mm error; conflicting anchors 2 m apart: splits the difference and stops, 0.02 m/s residual |
| overlapping pair with `AVBDIgnoreCollision3D` | 0 contact points, separation unchanged at 0.500 m; the same pair without it: 2 contact points and 1.019 m apart |
| breakable joint | `fracture_force` 20 breaks under a 40 rad/s spin (link released, keeps spinning at 10.5 rad/s); 2000 holds it (spin driven to 0). The break survives a rebuild |
| plank bridge: 22 joints + 11 ignores + 2 braces + 3 dropped loads | worst joint error 22 mm, 104 mm sag at mid-span, stays 5.66 m above the ground |
| 27-cell lattice on 9 springs (k = 600 N/m) | sag 1.03 mm vs 0.77 mm predicted by `m*g/(n*k)`; worst joint error 1.4 mm after 10 s |
| static driver teleported 4 m | rider dragged the full 3.9995 m, hitch re-converges in 437 steps, final error 0.57 mm, residual 0.3 mm/s |
| `substeps = 1` vs `2` on a joint + spring | same equilibrium within 0.94 mm |

## Limitations

* **Boxes only.** Collision is box-box SAT with face clipping, as in the reference.
  No convex hulls, meshes, capsules or spheres/triangles.
* **Broad phase is O(n²)** over bounding spheres, as upstream, and stays serial. Fine into
  the low thousands of bodies; a spatial hash or BVH is the obvious next step, and it is now
  the largest remaining serial cost.
* **Parallelism is limited by the constraint graph, not the core count.** Bodies coupled by
  a contact or a joint cannot be updated together, so a deep stack or a long chain has a
  narrow critical path (a stack is a path graph: 2 groups, half the bodies each). Scenes made
  of many weakly-coupled bodies scale much better than one tall pile.
* **Manifold creation is serial**, because it splices into the solver's and both bodies'
  force lists in an order that fixes the accumulation order of the solver. That ordering is
  what makes the results reproducible, so it is not parallelised.
* **No `PhysicsServer3D` backend.** These are ordinary scene nodes, not a drop-in
  replacement for Godot's physics server; Godot's own 3D physics can be disabled in the
  project settings (the demo does).
* Soft bodies are jointed-box lattices - there is no tetrahedral FEM/cloth solver.
* Static bodies move only by teleporting (their node transform); there is no kinematic
  velocity coupling, so a moving platform pushes through penetration correction only.
* Sleeping, CCD, and continuous contact caching across steps are not implemented.
* The extension targets desktop Linux (`linux.x86_64` in `demo/avbd.gdextension`);
  the xmake target itself is portable - `godotcpp4` builds for Windows/macOS/Android.

## Attribution

`src/avbd/` is a port of <https://github.com/savant117/avbd-demo3d> by Chris Giles
(MIT), with the reference's own license header retained in every file. Deliberate
deviations: the code is wrapped in `namespace avbd` and the rendering includes are
gone; `Solver::pick()` also tests static bodies so raycasts can hit the ground; the
API is extended with a `contactPointCount()` statistics hook; and `Rigid::constrainedTo()`
walks the body's own constraint list (`nextA`/`nextB`) instead of the solver-wide chain
(`next`). That last one is a performance fix, not a behaviour change - the two traversals
find the same constraints, but the upstream version costs O(all forces) for every broad-phase
pair test, which on the 528-box benchmark was 78% of the step time.

The core has since been brought up to modern C++ without changing a single number: the
arithmetic is still the arithmetic. The physics update was re-expressed as named phases
(`broadPhase`, `warmstartForces`, `colourGraph`, `warmstartBodies`, `solveIterations`,
`finishVelocities`) that mirror the paper; the six out-parameters of the primal update
became one `Block` (the 6x6 system a body solves, which now also knows how to solve and
apply itself); `defaultParams()` became member initialisers; the `PENALTY_*`/`COLLISION_*`
macros became `inline constexpr`; `using namespace std` is gone from the public header, as
are the C-style casts that aliased a vector's members as an array; and the intrusive lists
now have explicitly deleted copy/move operations, because copying one would share the list
heads. `std::exclusive_scan` does the colour-grouping prefix sum.

Two things that look like cleanups but are load-bearing, and were deliberately left alone:
the per-step colouring still rebuilds from scratch (it costs 0.01 ms for 500 bodies and its
first-fit assignment *is* the update order, hence part of the reproducibility guarantee),
and `Block` divides by `dt*dt` rather than multiplying by a reciprocal, because those are
not the same floating-point operation. Verified by the two recorded digests being
bit-identical before and after.

```bibtex
@article{Giles2025,
  author = {Chris Giles and Elie Diaz and Cem Yuksel},
  title = {Augmented Vertex Block Descent},
  journal = {ACM Transactions on Graphics (Proceedings of SIGGRAPH 2025)},
  year = {2025}, volume = {44}, number = {4},
}
```

Licensed under MIT; see `LICENSE`.

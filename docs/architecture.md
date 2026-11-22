# Source map

This is a reading guide, not a design essay. The invariants that matter most are
at the bottom.

## Solver core (`src/avbd/`)

| File | What lives there |
|---|---|
| `types.hpp` | `ShapeType`, `Shape`, `Rigid`, `Block`, and the shared constants. No solver state. |
| `manifold.hpp` | `Manifold`, the SoA `ManifoldStore`, `ForceRef`/`ForceAdjacency`, and the `collideShapes` declaration. |
| `forces.hpp` | `Force` and its concrete kinds: `Joint`, `GenericJoint`, `Spring`, `IgnoreCollision`. |
| `solver.h` | The `Solver` class plus the includes above; kept as the single umbrella include so existing callers need no changes. |
| `solver_internal.hpp` | File-local helpers shared by solver translation units (adjacency lookup, axis-lock masking). Not part of the public API. |
| `solver.cpp` | Construction and destruction only. |
| `solver_query.cpp` | `pick`, `collectAABB` — read-only queries over the last built broad-phase state. |
| `solver_step.cpp` | `step`, `stepInterleaved`, `runSubstep`, `rebaseContacts`, thread-pool setup — the sequencing of a step. |
| `solver_projection.cpp` | `projectContacts` — the second, Jacobi position-level solve. |
| `solver_broadphase.cpp` | LBVH build, pair filtering, sleep/wake on AABB overlap. |
| `solver_solve.cpp` | Warm start, the coloured Gauss-Seidel primal, dual updates, finish/update-set/colouring bookkeeping. |
| `manifold.cpp` | Narrow-phase data assembly, primal/dual contact solves, warm-start caching. |
| `collide.cpp` | Shape-pair collision (box SAT, sampled cylinder paths, sphere cases). |
| `joint.cpp`, `generic_joint.cpp`, `spring.cpp`, `force.cpp`, `rigid.cpp` | Constraint and body implementations. |
| `maths.h`, `linalg.h` | Linear algebra layer; `linalg.h` is third-party and should not be edited. |
| `bvh_broadphase.*`, `broadphase_pair.hpp`, `radixsort.hpp` | Broad-phase data structures. |
| `job_pool.hpp`, `list_range.hpp`, `fixed_stack.hpp` | Parallel dispatch and small containers. |
| `query/query_shape.*` | Solver-free ray/point/sweep queries, used by the Godot space state. |

## Godot layer (`src/server/`, `src/nodes/`)

One `Solver` per space; `AVBDPhysicsServer3D` owns the RID maps and rebuilds a
space's solver from the bodies Godot sends. `AVBDDirectSpaceState3D` and
`AVBDDirectBodyState3D` translate engine queries into the solver's Z-up frame.
The soft-body node pair in `src/nodes/` is separate from the server path.

## Invariants to preserve

- **Bit-exact determinism.** `xmake run avbd_core_test` prints FNV digests of
  the full body state; cross-thread equality is asserted in-process. Any change
  to a float expression's order, a loop order, or the `-ffp-contract=off`
  build flag changes those digests and must be treated as a behaviour change.
  `tools/verify_refactor.sh` checks the digests against
  `tools/baseline/digests.txt`.
- **Collision feature order.** Contact sampling order feeds manifold contact
  order, which feeds solve order. Deduplication must not reorder candidates.
- **SoA manifold store.** `ManifoldStore` keeps contacts in parallel arrays at a
  fixed stride; converting it to per-object records changes cache behaviour and
  the hot loops.
- **No feature removal.** The README's capability and limitation lists are the
  contract; simplifying code must not silently narrow either.

## Verification

```sh
tools/verify_refactor.sh          # core + query suites and digest comparison
```

Godot-side suites (`xmake run avbd_query_tests`, `avbd_solver_tests`, etc.)
require a Godot 4 binary and are not covered by that script.

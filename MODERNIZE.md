# Modernize branch: avbd-mo status

`modernize` replaces the legacy `src/avbd` solver with the imported core in
`src/modern` (see `src/modern/UPSTREAM.md` for provenance and sync rules). The
legacy core stays only for the soft-body nodes until Phase 4.

This file is the honest state of the migration, not a wish list.

## Working

- Godot physics-server main path: box / sphere / cylinder / capsule (cylinder
  proxy) / convex (AABB proxy) shapes, compound bodies (every enabled shape on a
  body becomes a `collision::ShapeInstance`), contacts, gravity and per-body
  gravity scale, damping, constant/one-shot forces and impulses, teleports,
  collision layers/masks/exceptions, explicit sleep, axis locks, direct body
  state, contact monitoring, and write-back through the normal sync callbacks.
- Space queries: ray, point, shape intersect, cast/rest, collision exceptions,
  areas and swept body motion (`body_test_motion`). Ray/point tests are exact
  for box/sphere/cylinder; shape-shape queries use conservative world AABBs.
- `avbd_modern_test`: 50/50 structural selftests, solver digest and collision
  digest bit-identical to the pinned upstream baselines. Wired into
  `run_all_tests` and CI (`avbd_modern_tests`).

## Known gaps (upstream, by design of this milestone)

- **Friction**: the imported core emits tangent rows but never updates their
  dual, which would turn them into rigid anchors. The adapter therefore forces
  contact materials to zero friction; boxes slide. `test_solver`'s friction
  cases fail, and a large pyramid cannot stand without friction.
- **Joints**: no joint constraint in the core yet. Joint records are stored and
  getters replay, but nothing is simulated (`test_solver`'s hard/world joint
  cases fail; `test_constraints` is expected to fail).
- **Automatic sleep**: only explicit `BODY_STATE_SLEEPING` freezes a body.
- **Speculative contacts / CCD / position projection**: not ported.
- **Capsule and convex** are documented proxies, not exact geometry.

Current `test_solver.gd`: 38 checks, 6 failures (3 friction, 2 joints, 1 impulse
travel distance). `test_server.gd` 30/30 and `test_queries.gd` 27/27 pass.

## Coordinates / conventions

Both sides are Y-up and the modern cylinder axis is local +Y, so the server path
has no axis permutation. The soft-body nodes keep their own Z-up conversion.

## Next

1. Implement friction in avbd-mo (dual update per manifold, tangent clamp to
   `mu * |lambda_n|`), re-record its baselines, re-copy.
2. Port joints as core constraints (`AxisConstraint`/generic 6-DOF).
3. Island-based sleep in the frame/world layer.
4. Then remove `src/avbd` and the soft-body dependence on it.

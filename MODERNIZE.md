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
- **Coulomb friction** (`7107cff`), following `savant117/avbd-demo3d`:
  per-row adaptive penalty, `lambda` stored as the clamped force,
  manifold-level dual (`F = K*C + lambda`, normal `<= 0`, tangents radially
  clamped to `mu*|F_n|`), first-order stabilization `C_eff = C - 0.99*C0`,
  warm-start decay `lambda *= alpha*gamma`, `K = clamp(K*gamma)`, and per-row
  primal force clamps with a lagged cone. Two additions are specific to this
  architecture's block-Jacobi gather: a mass-split Hessian weight and the
  per-row force clamp.
- That work also fixed a pre-existing sign error in
  `ContactConstraint::gradient_angular` (B-side torque was inverted). Without
  friction it showed up as slow wrong rolling; with a real dual it pumped
  energy every frame.
- `avbd_modern_test`: 50/50 structural selftests, solver digest and collision
  digest bit-identical to the pinned upstream baselines. Wired into
  `run_all_tests` and CI (`avbd_modern_tests`).

## Known gaps (upstream)

- **Joints**: no joint constraint in the core yet. Joint records are stored and
  getters replay, but nothing is simulated (`test_solver`'s hard/world joint
  cases fail; `test_constraints` is expected to fail). `hard_joint` used to pass
  on the frictionless tree only by accident: the links free-fell in formation
  because contacts were forced to zero friction.
- **Automatic sleep**: only explicit `BODY_STATE_SLEEPING` freezes a body.
- **Speculative contacts / CCD / position projection**: not ported.
- **Capsule and convex** are documented proxies, not exact geometry.

Current `test_solver.gd`: 38 checks, **3 failures** — all joints. Friction
(static creep 0.009 m, dynamic travel 4.3 m), stack height/penetration/drift
and both impulse checks pass. `test_server.gd` 30/30 and `test_queries.gd`
27/27 pass. `test_pyramid` still cannot hold a 420-box pile (joints and sleep
are absent, and a friction-only penalty model is not enough for that depth).

## Coordinates / conventions

Both sides are Y-up and the modern cylinder axis is local +Y, so the server path
has no axis permutation. The soft-body nodes keep their own Z-up conversion.

## Next

1. Port joints as core constraints (`AxisConstraint`/generic 6-DOF).
2. Island-based sleep in the frame/world layer.
3. Then remove `src/avbd` and the soft-body dependence on it.

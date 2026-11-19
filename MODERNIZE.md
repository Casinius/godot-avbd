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
- **Speculative contacts + CCD** (`8b9fa9a`): narrowphase generates true-gap
  contacts within `speculative_distance` and a pre-solve pass clamps closing
  speed to `(gap+margin)/dt`; fast bodies additionally get a translation-only
  swept TOI with momentum shared by inverse mass. The core default is 0, so
  pinned digests stay bit-identical; the server maps
  `physics/avbd/speculative_margin` (0.02 in the demo). Documented limits: the
  speculative pre-pass applies a central impulse (no torque), and CCD ignores
  rotation.
- **Contact stability** (`0c73d9a`): impact penetration is capped at 2 cm for
  the normal stabilization offset (the live bracket stays unclamped so deep
  overlaps keep their restoring gradient), contacts release their warm-started
  force when the fresh gap exceeds margin + 2 mm, and sticky anchors heal after
  2 mm of drift. This took `test_pyramid` from 7/9 to 9/9.
- **Joints** (`0c4853b`): a joint is decomposed into per-axis rows (locked,
  limited, free; linear and angular) with stabilization, force bounds and a
  penalty cap. The server maps pin, hinge, slider and generic 6DOF onto them,
  including world anchors. Cone-twist is deliberately not mapped: a swing cone
  is not a per-axis interval, and a wrong approximation would be worse than no
  constraint.
- `avbd_modern_test`: 50/50 structural selftests, solver digest and collision
  digest bit-identical to the pinned upstream baselines. Wired into
  `run_all_tests` and CI (`avbd_modern_tests`).

## Known gaps (upstream)

- **Cone-twist joint**: records and getters work, but no constraint is emitted
  (see above).
- **Automatic sleep**: only explicit `BODY_STATE_SLEEPING` freezes a body.
- **Position projection**: not ported. AVBD's position solve covers most of it;
  the remaining fidelity gap is the 40 s deep pile (below).
- **Deep-pile long run**: `test_pyramid` passes its 10 s gate (9/9), but a 40 s
  run still develops a late runaway. Deep-pile penetration in the 20–30 s
  window is ~14 cm, above the legacy-quality target.
- **Capsule and convex** are documented proxies, not exact geometry.

Current suites: `test_solver.gd` **38/38**, `test_constraints.gd` **11/11**,
`test_server.gd` **30/30**, `test_queries.gd` **27/27**, `test_load.gd` gate
PASS, `test_pyramid.gd` **9/9**.

## Coordinates / conventions

Both sides are Y-up and the modern cylinder axis is local +Y, so the server path
has no axis permutation. The soft-body nodes keep their own Z-up conversion.

## Next

1. Cone-twist joint (or a real cone constraint primitive).
2. Island-based sleep in the frame/world layer.
3. Revisit the 40 s deep-pile run with the new CCD path.
4. Then remove `src/avbd` and the soft-body dependence on it.

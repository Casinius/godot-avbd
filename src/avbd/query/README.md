# `avbd::query` - space-query extension layer

This folder is deliberately **outside** the teaching core. `solver.h`, `solver.cpp`,
`collide.cpp`, `maths.h` and the other solver files are unchanged; this module only reads
the public `avbd::Shape` description and calls the existing `collideShapes()` narrow phase.

## Why it exists

The solver itself only needs one narrow-phase entry (`collideShapes`) and one ray helper
(`Solver::pick`). Godot's space queries need a little more:

* exact ray/shape tests with a `t`, hit point, normal and "started inside" flag;
* a tight world AABB for a query shape;
* closest-point / containment for point and editor queries;
* a **swept** cast (`cast_motion`, `ShapeCast3D`, `CharacterBody3D` motion) that cannot
  tunnel through a target just because its final pose is past it;
* a penetration depth along the contact normal, so `rest_info` can rank candidates.

Putting those in a separate header keeps the solver readable and makes the added
convenience code easy to test on its own (`xmake run avbd_query_test`).

## Proxy rule

Queries use the **same proxy the solver simulates**:

| Godot shape | Query shape |
|---|---|
| Box | box |
| Sphere | sphere |
| Cylinder | Z-axis cylinder |
| Capsule | cylinder over the capsule radius and full height |
| Convex polygon | the points' local AABB as a box |
| Concave / heightmap / world boundary / separation ray | no query proxy (skipped) |

A raycast therefore cannot report a surface the solver does not collide against. When the
solver gains exact capsule/convex geometry, this module is the place to extend, without
touching the core.

## API

`query_shape.hpp`:

* `intersectRay(shape, ray, hit)` - box / sphere / cylinder, including caps and inside hits.
* `closestPoint(shape, point)` / `containsPoint(shape, point, margin)`.
* `shapeAabb(shape, min, max)`.
* `shapeContacts(a, b, contacts, basis, margin)` - the existing narrow phase.
* `sweepShape(base, motion, target, margin)` - bounded by the swept-AABB overlap interval,
  then sampled/bisected for a stable `safe` / `unsafe` pair.
* `supportRadius(shape, axis)` / `penetrationDepth(a, b, normal)` - rest-query ranking.

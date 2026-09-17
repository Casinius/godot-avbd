# AVBD for Godot 4

[AVBD](https://github.com/savant117/avbd-demo3d) (Augmented Vertex Block Descent) as a
Godot 4.7 physics engine, exposed as a GDExtension. Scenes built from Godot's own nodes
— `RigidBody3D`, `StaticBody3D`, `CollisionShape3D`, `Generic6DOFJoint3D` — are simulated
by the AVBD solver when the project selects it:

```ini
[physics]
3d/physics_engine="AVBD"
```

## Architecture: server first

The main path is the **PhysicsServer3D extension** (`src/server/`), the engine's own
extension point: the editor, scene files, `SceneTree` physics, and third-party addons
all keep working. Key properties:

- **Follows the game.** The solver steps in the engine's physics loop only. Pausing the
  `SceneTree` calls `PhysicsServer3D.set_active(false)` and the simulation freezes with
  everything else.
- **Standard nodes end to end.** Body state flows back through the Godot 4.x sync
  contract: the server calls the engine's callback with a `PhysicsDirectBodyState3D`
  (implemented by `AVBDDirectBodyState3D`), which is what drives `RigidBody3D`.
- **Space queries** (`World3D.direct_space_state`: `intersect_ray`, `intersect_point`,
  `intersect_shape`, `cast_motion`, `collide_shape`, `rest_info`, `RayCast3D`,
  `ShapeCast3D`) are answered by `AVBDDirectSpaceState3D`, using the same
  `collideShapes()` path the solver's broad phase uses.
- **Areas** run enter/exit monitoring each step and fire `Area3D`'s
  `body_entered`/`body_exited` callbacks with the engine's argument contract.
- **Solver tuning** lives in project settings under `physics/avbd/*`:
  `threads` (0 = auto), `iterations` (10), `alpha` (0.99), `beta_linear` (10000),
  `beta_angular` (100), `gamma` (0.999), `substeps` (1).

The only node classes left are the **soft-body pair** (`src/nodes/`): `AVBDSoftBody3D`
(a lattice of jointed boxes) and `AVBDSoftWorld3D`, which owns the solver soft bodies
need because they have no `PhysicsServer3D` representation yet. The old node layer
(`AVBDWorld3D`, `AVBDRigidBody3D`, `AVBDJoint3D`, `AVBDGeneric6DOFJoint3D`,
`AVBDSpring3D`, `AVBDIgnoreCollision3D`, `AVBDConstraint3D`) has been removed.

## Support matrix

| Feature | Status | Notes |
| --- | --- | --- |
| Rigid bodies (box/sphere/cylinder) | yes | one shape per body; extras are ignored with a warning |
| Collision layers/masks | yes | Godot pair semantics (`A.layer & B.mask`) |
| Collision exceptions | yes | `add_collision_exception_with` |
| Axis locks (linear/angular) | yes | locked axes resist gravity and impulses; joint coupling can leak within a step (approximation) |
| Sleeping (manual) | yes | `RigidBody3D.sleeping`; automatic idle-sleep is not implemented (the solver has no sleep trigger yet) |
| Direct space state queries | yes | ray/point/shape/cast/collide/rest_info |
| Pin / Hinge / Slider / ConeTwist joints | yes | mapped onto the solver's `GenericJoint`; off-axis hinge/slider/twist axes snap to the nearest principal axis |
| Generic6DOFJoint3D (limits + per-axis springs) | yes | Godot spells "locked" as limits 0..0 |
| Area monitoring | yes | enter/exit callbacks; area `gravity_override` etc. are recorded, not applied |
| Contact reporting | yes | `max_contacts_reported` caps what the server records |
| Scripted integration | yes | `_integrate_forces` via the force-integration callback |
| Damping (`linear_damp`, `angular_damp`) | recorded, not simulated | the solver has no damping term; getters replay the values |
| Automatic sleeping | no | manual `sleeping` only |
| `_body_test_motion` (`CharacterBody3D`) | no | motion queries are not implemented; `CharacterBody3D` will not move |
| Soft-body server API | no | soft bodies run through the `AVBDSoftWorld3D` node pair instead |
| Capsule / convex / concave / heightmap shapes | no | only box, sphere, cylinder |
| CCD | recorded, not simulated | fast tunnelling is mitigated only by the solver's contact margin |
| Area gravity overrides / Ray picking flags | recorded | not consulted by queries this round |

## Building

```bash
xmake f -m release && xmake build avbd        # GDExtension (copied to demo/bin/)
xmake build avbd_core_test                    # standalone solver test binary

# Import the demo project once so the extension is registered
godot --headless --path demo --import
```

Sanitizers: `xmake f -m debug --asan=y --ubsan=y && xmake` (each mode needs its own
build directory).

## Testing

```bash
xmake run avbd_core_test        # Tier 1: C++ solver tests (no engine; 83 checks)
xmake run check_api             # Tier 2: server API compliance audit
xmake run avbd_solver_tests     # Tier 3: solver behaviour in-engine (test_solver.gd)
xmake run avbd_constraint_tests # Tier 4: constraint composition in-engine (test_constraints.gd)
xmake run avbd_load_tests       # Load gate: class registration + 200-box pyramid
xmake run run_all_tests         # Everything above, in order
```

Every Godot-side suite prints `physics_engine=AVBD` (or fails the engine check) to prove
it is exercising the server path. The suites live in `demo/tests/` and share
`avbd_harness.gd` (scenario declaration, assertions, `--digest` reproducibility mode);
each scenario builds its scene from standard nodes and steps the real physics loop.
CI runs `xmake run run_all_tests`; see `.github/workflows/physics_tests.yml`.

## Demos

`demo/scenes/` — all standard nodes, driven by the AVBD server:

- `stack.tscn` — a 10-box stack (default scene)
- `friction.tscn` — two boxes on a 20-degree ramp (`PhysicsMaterial.friction` 0.5 vs 0.05)
- `pyramid.tscn` — a 210-box pyramid
- `bridge.tscn` — planks hinged by `Generic6DOFJoint3D` between two posts
- `soft_body.tscn` — three soft lattices via `AVBDSoftWorld3D`

Controls (windowed): left click punches along the camera ray (a direct-space-state
`intersect_ray`), **Space** pauses the `SceneTree` (the server freezes with it),
**R** resets every body to its scene pose.

## License

MIT (see LICENSE); the solver core carries the original authors' notice.

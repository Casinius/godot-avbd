# avbd-mo snapshot

This tree is a snapshot of `/home/cyan/avbd-mo/src`, imported to replace the legacy
`src/avbd` solver on the `modernize` branch.

- Source HEAD: `a3d99303ff6996446f71e9533bb799ca80e09aa2` (`manifold`)
- Source working tree: dirty. The staged collision-layer work was included:
  `collision_digest.hpp`, `convex_hull.*`, `convex_mesh.*`, `mass_properties.hpp`,
  `shape_instance.hpp`, `shape_table.hpp`, plus the modified broadphase/narrowphase/
  manifold/collision_world/shape/main files.
- Snapshot tree SHA-256: `31a80374eb9afe3cfe3bfafa345b1e02ee9ff90034947141823905f46233590d`
- Baseline data copied to `src/modern/tests/`:
  `baseline.txt`, `collision_observable.txt` (the solver/collision digests this snapshot
  was verified against).

## Rules

- Keep this directory as close to upstream as possible. Godot-specific code belongs in
  `src/godot`, `src/server`, or `src/nodes`.
- If a change is genuinely generic (collision filtering, friction, joints, sleep),
  implement and test it in `/home/cyan/avbd-mo` first, then re-copy the snapshot and
  update the hashes above.

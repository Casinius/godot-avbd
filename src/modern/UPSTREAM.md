# avbd-mo snapshot

This tree is a snapshot of `/home/cyan/avbd-mo/src`, imported to replace the legacy
`src/avbd` solver on the `modernize` branch.

- Source HEAD: `4450b74` (`manifold`) — "add collision layer: compound shapes, convex
  meshes, digests, pair filter". The earlier snapshot was taken from `a3d9930` plus its
  staged collision work; that work is now committed upstream, so this tree corresponds
  to a clean commit.
- Snapshot tree SHA-256 (everything under `src/modern` except this file):
  `99d1ecb823736ee01ebe1722ef9af1c52c8d2b8b0f996276fbfed333ab77f711`
- Baseline data copied to `src/modern/tests/`:
  `baseline.txt`, `collision_observable.txt` (the solver/collision digests this snapshot
  was verified against).
- The snapshot builds with `-Wshadow` warnings at `phys_core/state.hpp:518` (an `int r`
  shadowing a parameter). Left as-is: it is upstream's code and the warning is harmless.

## Rules

- Keep this directory as close to upstream as possible. Godot-specific code belongs in
  `src/godot`, `src/server`, or `src/nodes`.
- If a change is genuinely generic (collision filtering, friction, joints, sleep),
  implement and test it in `/home/cyan/avbd-mo` first, then re-copy the snapshot and
  update the hashes above.

# avbd-mo snapshot

This tree is a snapshot of `/home/cyan/avbd-mo/src`, imported to replace the legacy
`src/avbd` solver on the `modernize` branch.

- Source HEAD: `04d2833` (`manifold`) — "fix block assembly symmetry; configurable
  AVBD k0; expose convex mesh store". Earlier snapshots were `4450b74` (pair filter)
  and `a3d9930` plus staged collision work; this tree corresponds to a clean commit.
- Snapshot tree SHA-256 (everything under `src/modern` except this file):
  `b5f2d64e029a36e5904be8a0017d6929575a751b432f77ec026deceee7ce77f9`
- Baseline data copied to `src/modern/tests/`:
  `baseline.txt`, `collision_observable.txt` (the solver/collision digests this snapshot
  was verified against).
- The snapshot builds with `-Wshadow` warnings in `phys_core/state.hpp` (an `int r`
  shadowing a parameter in `solve_body`). Left as-is: it is upstream's code and the
  warning is harmless.
- `phys_core/policy.hpp` now exposes `AVBDStep::k0`; the Godot adapter maps its
  `physics/avbd/beta_linear` setting onto it. Defaults are unchanged, so the pinned
  baselines still match.

## Rules

- Keep this directory as close to upstream as possible. Godot-specific code belongs in
  `src/godot`, `src/server`, or `src/nodes`.
- If a change is genuinely generic (collision filtering, friction, joints, sleep),
  implement and test it in `/home/cyan/avbd-mo` first, then re-copy the snapshot and
  update the hashes above.

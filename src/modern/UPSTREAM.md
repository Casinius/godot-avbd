# avbd-mo snapshot

This tree is a snapshot of `/home/cyan/avbd-mo/src`, imported to replace the legacy
`src/avbd` solver on the `modernize` branch.

- Source HEAD: `7107cff` (`manifold`) — "implement Coulomb friction and fix contact
  angular gradient sign". Earlier snapshots were `04d2833` (k0 / block assembly),
  `4450b74` (pair filter) and `a3d9930` plus staged collision work; this tree
  corresponds to a clean commit.
- Snapshot tree SHA-256 (everything under `src/modern` except this file):
  `712fecaeeaeef1e7dbe74619009dfe6ad41e9b3f456ae9b54c3c5a172ef3b93f`
  Recipe: `find src/modern -type f ! -name UPSTREAM.md -print0 | sort -z |
  xargs -0 sha256sum | sha256sum`.
- Baseline data copied to `src/modern/tests/`:
  `baseline.txt`, `collision_observable.txt` (the solver/collision digests this snapshot
  was verified against). The solver baseline was re-recorded when contact duals became
  real; the collision digest is unchanged.
- The snapshot builds with `-Wshadow` warnings in `phys_core/state.hpp` (an `int r`
  shadowing a parameter in `solve_body`). Left as-is: it is upstream's code and the
  warning is harmless.
- `phys_core/policy.hpp` exposes `AVBDStep::k0`; the Godot adapter maps its
  `physics/avbd/beta_linear` setting onto it. `k0` doubles as the contact penalty
  ceiling for the adaptive ramp.

## Rules

- Keep this directory as close to upstream as possible. Godot-specific code belongs in
  `src/godot`, `src/server`, or `src/nodes`.
- If a change is genuinely generic (collision filtering, friction, joints, sleep),
  implement and test it in `/home/cyan/avbd-mo` first, then re-copy the snapshot and
  update the hashes above.

# avbd-mo snapshot

This tree is a snapshot of `/home/cyan/avbd-mo/src`, imported to replace the legacy
`src/avbd` solver on the `modernize` branch.

- Source HEAD: `0c73d9a` (`manifold`) — "contact stability: penetration cap,
  release, stick-anchor heal". Earlier snapshots were `0c4853b` (joints),
  `7107cff` (friction + gradient sign), `04d2833` (k0 / block assembly),
  `4450b74` (pair filter) and `a3d9930` plus staged collision work.
- The contact-stability part of this snapshot:
  * `phys_core/constraint.hpp` — normal-row penetration cap
    `kContactMaxPenetration = 0.02f` (legacy `CONTACT_MAX_PENETRATION`), applied to
    the step-start stabilization offset on impact steps; `initial_c_raw` keeps the
    unclamped C0 so `evaluate()` retains the live restoring bracket (the literal
    "clamp the live C" form makes a deeply overlapped pair lose all support and sink).
  * `collision/manifold.{hpp,cpp}` — contact release (`kContactReleaseGap = 2 mm`)
    and stick-anchor heal (`kStickDriftHeal = 2 mm`), matching legacy
    `CONTACT_RELEASE_GAP` / `STICK_DRIFT_HEAL`; `ContactPatch::margin` carries the
    emit margin into `ManifoldStore::match`.
  * No pinned solver/collision digest changes: the cap is inert on the pinned
    scenes, and the release path never fires in the modern (non-speculative)
    narrowphase.
- Snapshot tree SHA-256 (everything under `src/modern` except this file):
  `cf9951f08e783a19d03640f15f3ce9fabc893a5d9d4affa8ac5f1f44c3f53e93`
  Recipe: `find src/modern -type f ! -name UPSTREAM.md -print0 | sort -z |
  xargs -0 sha256sum | sha256sum`.
- Baseline data copied to `src/modern/tests/`:
  `baseline.txt`, `collision_observable.txt` (the solver/collision digests this snapshot
  was verified against). The solver baseline was re-recorded when contact duals became
  real; adding joints did not change it (the pinned scenes are joint-free), and the
  collision digest is unchanged.
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

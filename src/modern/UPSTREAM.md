# avbd-mo snapshot

This tree is a snapshot of `/home/cyan/avbd-mo/src`, imported to replace the legacy
`src/avbd` solver on the `modernize` branch.

- Source HEAD: `4cd21cf` (`manifold`) — "BPXPBD: opt-in velocity-neutral
  end-of-step projection" (on top of `feca2f9` cap/friction, `824b37b`
  contacts and `e6f995a` skeleton). Earlier snapshots were `8b9fa9a` (speculative
  contacts/CCD), `0c73d9a` (contact stability), `0c4853b` (joints), `7107cff`
  (friction + gradient sign), `04d2833` (k0 / block assembly), `4450b74`
  (pair filter) and `a3d9930` plus staged collision work.
- The speculative/CCD part of this snapshot (default `speculative_distance = 0`,
  so the pinned solver/collision baselines are still bit-identical):
  * `collision/contact_constraints.hpp` — `ContactEmitOptions::speculative_distance`
    is a generation distance, distinct from the resting `margin`. Contacts with
    `depth < 0` (a true gap) are not emitted as solver rows: the normal row's
    `K*C` would attract while `C = gap + margin > 0`.
  * `collision/narrowphase.{hpp,cpp}` — the exact overlap path is unchanged; a
    second pass generates `depth = -gap` contacts at the closest features within
    `speculative_distance` (sphere/sphere, sphere/polytope, polytope face/edge/
    support). Face/edge SAT now has a gap-tolerant variant; `sweep_shapes_toi`
    provides a conservative translation-only TOI for CCD (rotation ignored).
  * `collision/collision_world.{hpp,cpp}` — broadphase proxy AABBs expand by the
    generation distance (plus a swept proxy for fast bodies); `clamp_fast_bodies`
    is the momentum-conserving swept fallback (earliest TOI, relative-velocity
    scaling, small skin), and `apply_speculative_velocity_clamp` reduces closing
    speed to `(gap + margin)/dt` before the position solve using the full
    effective mass (inverse mass + world inverse inertia). The speculative
    pre-pass applies a **central** impulse; torque is left to exact contact rows
    because projected speculative anchors otherwise twist stacks apart.
  * `collision/aabb.hpp` — `swept_aabb_fraction` (conservative swept-AABB helper,
    rotation ignored).
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
    scenes, and the release path only fires for speculative gaps (which are off
    by default). The snapshot still matches `tests/baseline.txt` and
    `tests/collision_observable.txt` bit-for-bit.
- Snapshot tree SHA-256 (everything under `src/modern` except this file):
  `3a0efa7a768934de267e4ecacd7413bd922c30b4c679d99a6f96062ed263a228`
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

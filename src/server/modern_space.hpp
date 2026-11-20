/*
 * ModernSpace - the Godot-facing adapter around the imported avbd-mo core.
 *
 * Part of godot-avbd, on the `modernize` branch.
 *
 * What this owns
 * --------------
 * One instance per Godot space. It holds the whole modern pipeline:
 *
 *     State (pose/velocity/inertia/external) + ConstraintList
 *     collision::CollisionWorld (shape table -> broadphase -> narrowphase -> manifolds)
 *     std::variant<AVBDStep, BPXPBDStep> + SolveSchedule + a task::ThreadPoolExecutor
 *
 * and nothing Godot-specific. The server translates its records into ModernBodyDesc
 * once per rebuild, then drives this object with begin_tick/step.
 *
 * `SolverKind` selects the policy at runtime (project setting `physics/avbd/solver`).
 * AVBD remains the default; switching zeroes all dual state and clears the collision
 * world's warm start, because AVBD's lambda is a force while B-PXPBD's is a position
 * multiplier (different units, not interchangeable).
 *
 * Design notes that the rest of the server depends on
 * ---------------------------------------------------
 * - **Coordinates are Godot's.** The modern core's cylinder axis is +Y and its gravity
 *   is just a default uniform acceleration, so no Z-up/Y-up permutation is needed
 *   anywhere on this path. (The legacy soft-body nodes still run their own Z-up core.)
 * - **Body indices are stable.** `rebuild` sorts by server body id, so the state/solver
 *   order - and therefore the Gauss-Seidel order and every digest - is a function of
 *   the scene, not of hash-map layout.
 * - **Pose is the centre of mass.** Godot body origins are not necessarily the COM of
 *   their compound shapes, so rebuild recentres the shape instances and keeps the
 *   body-origin offset (`com_local`). set_body_transform/body_origin convert, and the
 *   server never has to think about it.
 * - **Layer/mask/exceptions are enforced by collision::CollisionWorld::set_pair_filter**,
 *   not by a second broad phase. The filter reads per-index layer/mask tables and a
 *   per-index sorted exception list, so the step path allocates nothing.
 * - **Friction is not implemented by the imported core yet.** Tangent contacts are
 *   emitted but their dual is inert; boxes will slide. Documented, not hidden.
 */

#ifndef AVBD_MODERN_SPACE_HPP
#define AVBD_MODERN_SPACE_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <variant>
#include <vector>

#include "collision/collision_world.hpp"
#include "collision/contact_constraints.hpp"
#include "collision/mass_properties.hpp"
#include "collision/shape_instance.hpp"
#include "phys_core/joint.hpp"
#include "phys_core/policy.hpp"
#include "phys_core/solver_phase.hpp"
#include "phys_core/state.hpp"
#include "task.hpp"

namespace godot {

// One body, as the server describes it to the modern core. Pure modern types: this
// header must stay free of godot-cpp so the adapter can be reasoned about (and tested)
// on its own.
struct ModernBodyDesc {
    std::uint64_t id = 0; // server body id, used only for diagnostics/exceptions
    bool is_static = false;
    double mass = 1.0; // kg; ignored when is_static
    Vec3 position{0}; // body origin, world space
    Quat rotation = identity_quat();
    Vec3 linear_velocity{0};
    Vec3 angular_velocity{0};
    std::uint32_t layer = 1;
    std::uint32_t mask = 1;
    // Godot BodyAxis bit order: linear x/y/z = 1/2/4, angular = 8/16/32.
    std::uint32_t axis_locks = 0;
    std::vector<collision::ShapeInstance> shapes;
    // Server ids this body must never collide with (Godot's exception list).
    std::vector<std::uint64_t> exceptions;
    // Convex meshes are registered by the server before rebuild; index 0 is reserved
    // for "no mesh". Kept out of ShapeInstance so the desc stays copyable/cheap.
};

// One joint axis, as the server describes it. Axis i is the i-th axis of the joint
// frame (frame_a / frame_b rotations, usually the identity); Free axes emit no row.
struct ModernJointAxis {
    JointAxisMode mode = JointAxisMode::Free;
    Scalar_t lower = 0; // metres for linear, radians for angular
    Scalar_t upper = 0;
    Scalar_t target = 0;
};

// One joint, as the server describes it to the modern core. Body indices are ModernSpace
// state indices; body_a < 0 means the joint's A side is the world, which is represented
// by the persistent static world body (see world_body_index()).
struct ModernJointDesc {
    std::uint64_t id = 0; // server joint id, used only to keep rebuild order deterministic
    std::int32_t body_a = -1;
    std::int32_t body_b = -1;
    Vec3 anchor_a{0}, anchor_b{0}; // body-local anchors
    Quat frame_a = identity_quat();
    Quat frame_b = identity_quat();
    std::array<ModernJointAxis, 3> linear{};
    std::array<ModernJointAxis, 3> angular{};
};

// Policy selection. AVBD is the default; B-PXPBD is opt-in via the project setting.
enum class SolverKind {
    AVBD = 0,
    BPXPBD = 1,
};

class ModernSpace {
public:
    ModernSpace();
    ~ModernSpace();

    ModernSpace(const ModernSpace &) = delete;
    ModernSpace &operator=(const ModernSpace &) = delete;

    // Replace the whole body set. Sorts by desc.id, builds shapes, mass properties,
    // collision filter tables and the collision world. Destroys contact warm start.
    // Also reserves one trailing static body for joints whose A side is the world, so a
    // joint always has two valid state indices.
    void rebuild(std::vector<ModernBodyDesc> bodies);
    void clear();

    // Replace every joint row. Rows are kept in `constraints_` in front of the contact
    // rows (collision_.sync only erases/re-emits contacts at the tail), so repeated steps
    // do not rebuild them and keep their warm-started dual state. Bumps topology_version_
    // only when the row identities change. Dual state is carried over for rows that keep
    // the same (bodies, kind, axis, mode).
    void rebuild_joints(std::vector<ModernJointDesc> joints);

    // Index of the persistent static world body. Valid after rebuild(); it is the state
    // index just past the last real body and is never part of ids_/collision pairs.
    [[nodiscard]] std::size_t world_body_index() const noexcept { return world_index_; }

    [[nodiscard]] std::size_t size() const noexcept { return ids_.size(); }
    [[nodiscard]] bool empty() const noexcept { return ids_.empty(); }
    [[nodiscard]] std::uint64_t body_id(std::size_t index) const noexcept { return ids_[index]; }
    // -1 when the id is not part of this space.
    [[nodiscard]] std::int32_t index_of(std::uint64_t body_id) const noexcept;

    // --- engine -> solver writes -------------------------------------------------
    // Positions are **body origins** in Godot space; the COM offset is applied here.
    void set_body_transform(std::size_t index, Vec3 origin, Quat rotation);
    void set_linear_velocity(std::size_t index, Vec3 velocity);
    void set_angular_velocity(std::size_t index, Vec3 omega);
    // Instant velocity changes, in world space. `world_offset` is the contact point
    // relative to the body origin (Godot semantics for apply_impulse).
    void apply_central_impulse(std::size_t index, Vec3 impulse);
    void apply_impulse(std::size_t index, Vec3 impulse, Vec3 world_offset);
    void apply_torque_impulse(std::size_t index, Vec3 impulse);
    void set_gravity_scale(std::size_t index, double scale);
    void set_damping(std::size_t index, double linear, double angular);
    // Explicit sleep (Godot's BODY_STATE_SLEEPING). A sleeping body keeps its pose and
    // zero velocity; it still participates in contacts as a rigid obstacle. Automatic
    // idle-sleep does not exist in the modern core yet.
    void set_sleeping(std::size_t index, bool sleeping);

    // --- solver -> engine reads --------------------------------------------------
    // Body origin (not COM), so the server can write Godot transforms directly.
    [[nodiscard]] Vec3 body_origin(std::size_t index) const;
    [[nodiscard]] Quat body_rotation(std::size_t index) const;
    [[nodiscard]] Vec3 linear_velocity(std::size_t index) const;
    [[nodiscard]] Vec3 angular_velocity(std::size_t index) const;
    [[nodiscard]] Vec3 com_local(std::size_t index) const noexcept { return com_local_[index]; }
    [[nodiscard]] bool is_static(std::size_t index) const noexcept { return state_.inertial.is_static(index); }

    // --- one tick ----------------------------------------------------------------
    // begin_tick resets the external accumulator to `gravity` and stores dt. Per-body
    // forces/torques go through add_force/add_torque afterwards; step() then runs the
    // full core pipeline. Splitting these is deliberate: frame::step_frame resets the
    // accumulator inside itself and cannot accept injected forces.
    void begin_tick(Vec3 gravity, Scalar_t dt);
    void add_force(std::size_t index, Vec3 force);
    void add_torque(std::size_t index, Vec3 torque);
    void step(int iterations, int threads);

    // --- tuning ------------------------------------------------------------------
    // Select the solver policy. No-op when the kind is unchanged. On a change it
    // zeroes every constraint's dual state (lambda/initial/initial_c/penalty where
    // applicable) and clears the collision world's warm start: AVBD's lambda is a
    // force, B-PXPBD's is a position multiplier, so they must not cross-inject.
    void set_solver(SolverKind kind);
    [[nodiscard]] SolverKind solver_kind() const noexcept;

    // Type-aware stiffness. AVBD maps (k0, beta) onto its penalty ramp; B-PXPBD uses
    // `zeta` for the blend term and `contact_compliance` (1/K, m/N) for contact rows.
    void set_policy_stiffness(Scalar_t k0, Scalar_t beta, Scalar_t zeta = Scalar_t(0.5),
            Scalar_t contact_compliance = Scalar_t(0)) noexcept;
    // B-PXPBD only: velocity-neutral post-solve position projection. AVBD ignores
    // this (its implicit contact solve already handles depenetration), so passing
    // values while AVBD is selected is a no-op. `enabled == false` or
    // `iterations <= 0` turns the pass off. Must be called after set_solver(),
    // which constructs a fresh BPXPBDStep on a switch.
    void set_contact_projection(bool enabled, int iterations, Scalar_t beta, Scalar_t slop,
            Scalar_t max_correction) noexcept;
    // Resting overlap the normal constraints target (metres). Smaller = less visible
    // penetration at the cost of a stiffer contact solve.
    void set_contact_margin(Scalar_t margin) noexcept;
    // Contact generation distance (metres). >0 turns on speculative contacts and the
    // conservative swept-CCD fallback. 0 keeps the core's exact-overlap-only behaviour.
    void set_speculative_distance(Scalar_t distance) noexcept;
    [[nodiscard]] Scalar_t speculative_distance() const noexcept {
        return collision_.speculative_distance();
    }
    [[nodiscard]] int threads_used() const noexcept;

    // --- observations ------------------------------------------------------------
    [[nodiscard]] const State &state() const noexcept { return state_; }
    [[nodiscard]] std::span<const collision::Manifold> manifolds() const noexcept {
        return collision_.manifolds();
    }
    [[nodiscard]] const collision::CollisionWorldStats &collision_stats() const noexcept {
        return collision_.stats();
    }
    [[nodiscard]] const collision::CollisionWorld &collision_world() const noexcept {
        return collision_;
    }
    [[nodiscard]] double last_step_ms() const noexcept { return last_step_ms_; }
    [[nodiscard]] std::uint64_t topology_version() const noexcept { return topology_version_; }

    // Convex meshes: the server registers them before rebuild (convex support is not
    // wired into the Godot shape path yet, but the registry lives here so it can be).
    collision::ConvexMeshId add_convex_mesh(collision::ConvexMesh mesh);

private:
    static bool pair_filter_cb(void *user, std::size_t body_a, std::size_t body_b);
    static bool body_active_cb(void *user, std::size_t body);
    void ensure_executor(int threads);
    // B-PXPBD only: inelastic normal-velocity projection on the manifolds generated
    // at the predicted pose. Runs before the position solve, fixed order, no atomics.
    void project_contact_normal_velocity();
    void apply_mass_properties(std::size_t index, const ModernBodyDesc &desc,
            std::span<const collision::ShapeInstance> shapes);
    void rebuild_filter_tables(const std::vector<ModernBodyDesc> &bodies);

    State state_{};
    ConstraintList constraints_{};
    std::vector<Joint> joints_{};
    std::size_t world_index_ = 0;
    collision::CollisionWorld collision_{};
    // B-PXPBD is the primary solver now; AVBD stays selectable as a legacy policy.
    std::variant<AVBDStep, BPXPBDStep> policy_{BPXPBDStep{}};
    SolveSchedule schedule_{};
    std::unique_ptr<task::ThreadPoolExecutor> exec_{};
    int exec_threads_ = -1;

    std::vector<std::uint64_t> ids_;
    std::vector<Vec3> com_local_;
    std::vector<std::uint32_t> layers_;
    std::vector<std::uint32_t> masks_;
    std::vector<std::vector<std::uint32_t>> exceptions_;
    std::vector<double> gravity_scale_;
    std::vector<double> linear_damp_;
    std::vector<double> angular_damp_;
    std::vector<std::uint32_t> axis_locks_;
    std::vector<bool> sleeping_;

    std::uint64_t topology_version_ = 0;
    double last_step_ms_ = 0.0;
};

} // namespace godot

#endif // AVBD_MODERN_SPACE_HPP

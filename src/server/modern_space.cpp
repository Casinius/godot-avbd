#include "server/modern_space.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <glm/matrix.hpp>
#include <type_traits>
#include <unordered_map>

namespace godot {

namespace {

// Godot BodyAxis bit order -> local axis index.
constexpr std::uint32_t kLockLinearX = 1u << 0;
constexpr std::uint32_t kLockLinearY = 1u << 1;
constexpr std::uint32_t kLockLinearZ = 1u << 2;
constexpr std::uint32_t kLockAngularX = 1u << 3;
constexpr std::uint32_t kLockAngularY = 1u << 4;
constexpr std::uint32_t kLockAngularZ = 1u << 5;

void zero_axis(Mat3 &m, int axis) {
    for (int j = 0; j < 3; ++j) {
        m[axis][j] = Scalar_t(0);
        m[j][axis] = Scalar_t(0);
    }
}
} // namespace

ModernSpace::ModernSpace() = default;
ModernSpace::~ModernSpace() = default;

void ModernSpace::clear() {
    state_.clear();
    constraints_.clear();
    joints_.clear();
    world_index_ = 0;
    collision_.clear();
    ids_.clear();
    com_local_.clear();
    layers_.clear();
    masks_.clear();
    exceptions_.clear();
    gravity_scale_.clear();
    linear_damp_.clear();
    angular_damp_.clear();
    axis_locks_.clear();
    sleeping_.clear();
    topology_version_ = 0;
    last_step_ms_ = 0.0;
}

collision::ConvexMeshId ModernSpace::add_convex_mesh(collision::ConvexMesh mesh) {
    return collision_.add_convex_mesh(std::move(mesh));
}

std::int32_t ModernSpace::index_of(std::uint64_t body_id) const noexcept {
    const auto it = std::lower_bound(ids_.begin(), ids_.end(), body_id);
    if (it == ids_.end() || *it != body_id)
        return -1;
    return static_cast<std::int32_t>(it - ids_.begin());
}

void ModernSpace::rebuild_filter_tables(const std::vector<ModernBodyDesc> &bodies) {
    const std::size_t n = bodies.size();
    layers_.assign(n, 1);
    masks_.assign(n, 1);
    exceptions_.assign(n, {});
    gravity_scale_.assign(n, 1.0);
    linear_damp_.assign(n, 0.0);
    angular_damp_.assign(n, 0.0);
    axis_locks_.assign(n, 0);
    sleeping_.assign(n, false);

    std::unordered_map<std::uint64_t, std::uint32_t> index_of_id;
    index_of_id.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        layers_[i] = bodies[i].layer;
        masks_[i] = bodies[i].mask;
        index_of_id.emplace(bodies[i].id, static_cast<std::uint32_t>(i));
    }

    // Godot exceptions are per-object and effectively symmetric. Resolve both
    // directions here so the step-time filter only does two binary searches.
    for (std::size_t i = 0; i < n; ++i) {
        for (const std::uint64_t other_id : bodies[i].exceptions) {
            const auto it = index_of_id.find(other_id);
            if (it == index_of_id.end() || it->second == i)
                continue;
            const std::uint32_t j = it->second;
            exceptions_[i].push_back(j);
            exceptions_[j].push_back(static_cast<std::uint32_t>(i));
        }
    }
    for (auto &list : exceptions_) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }
}

bool ModernSpace::pair_filter_cb(void *user, std::size_t a, std::size_t b) {
    const auto *self = static_cast<const ModernSpace *>(user);
    if (a >= self->layers_.size() || b >= self->layers_.size())
        return true;
    // Godot semantics: collide when either object's layer is in the other's mask.
    if ((self->layers_[a] & self->masks_[b]) == 0 &&
        (self->layers_[b] & self->masks_[a]) == 0)
        return false;
    const auto &ea = self->exceptions_[a];
    if (std::binary_search(ea.begin(), ea.end(), static_cast<std::uint32_t>(b)))
        return false;
    const auto &eb = self->exceptions_[b];
    return !std::binary_search(eb.begin(), eb.end(), static_cast<std::uint32_t>(a));
}

bool ModernSpace::body_active_cb(void *user, std::size_t body) {
    const auto *self = static_cast<const ModernSpace *>(user);
    // The persistent world body (index == sleeping_.size()) is static and always
    // "active"; only explicit sleepers are treated as immovable by the CCD pass.
    if (body >= self->sleeping_.size())
        return true;
    return !self->sleeping_[body];
}

void ModernSpace::apply_mass_properties(std::size_t index, const ModernBodyDesc &desc,
        std::span<const collision::ShapeInstance> shapes) {
    const bool static_body = desc.is_static || !(desc.mass > 0.0);

    collision::MassProperties mp;
    const collision::MassStatus status =
            collision::body_mass_properties(shapes, collision_.convex_meshes(), Scalar_t(1), mp);

    if (status != collision::MassStatus::Ok || !(mp.mass > 0)) {
        // No usable geometry: keep the body in the simulation as a point mass rather
        // than freezing it. Its shape list is empty, so it contributes no contacts.
        const Scalar_t m = static_cast<Scalar_t>(std::max(desc.mass, 1e-3));
        com_local_[index] = Vec3{0};
        state_.inertial.set_inv_mass(index, static_body ? Scalar_t(0) : Scalar_t(1) / m);
        state_.inertial.set_local_inv_inertia(index,
                static_body ? Mat3(0) : collision::diag3(Vec3{Scalar_t(2.5) / m}));
        return;
    }

    com_local_[index] = mp.com;
    if (static_body) {
        state_.inertial.set_inv_mass(index, 0);
        state_.inertial.set_local_inv_inertia(index, Mat3(0));
        return;
    }

    // Shape properties are computed with density 1; rescale so the body's authored
    // mass is what the solver sees while the inertia *ratio* stays geometric.
    const Scalar_t scale = static_cast<Scalar_t>(desc.mass) / mp.mass;
    const Mat3 inertia = mp.inertia * scale;
    Mat3 inv_inertia = glm::inverse(inertia);

    if (desc.axis_locks & kLockAngularX)
        zero_axis(inv_inertia, 0);
    if (desc.axis_locks & kLockAngularY)
        zero_axis(inv_inertia, 1);
    if (desc.axis_locks & kLockAngularZ)
        zero_axis(inv_inertia, 2);

    // A fully linear-locked body cannot translate. The solver's mass is a scalar, so
    // partial linear locks are not representable; the server documents that limit.
    const std::uint32_t linear = desc.axis_locks & (kLockLinearX | kLockLinearY | kLockLinearZ);
    const bool all_linear_locked = linear == (kLockLinearX | kLockLinearY | kLockLinearZ);

    state_.inertial.set_inv_mass(index, all_linear_locked ? 0 : Scalar_t(1) / static_cast<Scalar_t>(desc.mass));
    state_.inertial.set_local_inv_inertia(index, inv_inertia);
}

void ModernSpace::rebuild(std::vector<ModernBodyDesc> bodies) {
    std::sort(bodies.begin(), bodies.end(),
            [](const ModernBodyDesc &a, const ModernBodyDesc &b) { return a.id < b.id; });

    clear();
    const std::size_t n = bodies.size();
    // One extra state slot is the static world body that world-anchored joints use as
    // their A side. It has no id, no shapes and never takes part in collision; it is
    // simply a fixed frame the constraint math can read like any other body.
    world_index_ = n;
    state_.resize(n + 1);
    ids_.resize(n);
    com_local_.resize(n, Vec3{0});
    collision_.resize(n + 1);
    rebuild_filter_tables(bodies);

    for (std::size_t i = 0; i < n; ++i) {
        const ModernBodyDesc &desc = bodies[i];
        ids_[i] = desc.id;
        gravity_scale_[i] = 1.0; // overwritten by the server through set_gravity_scale
        linear_damp_[i] = 0.0;
        angular_damp_[i] = 0.0;
        axis_locks_[i] = desc.axis_locks;
        sleeping_[i] = false;

        std::vector<collision::ShapeInstance> shapes = desc.shapes;
        apply_mass_properties(i, desc, shapes);

        // Recentre the shape instances around the COM; the pose is the COM from here on.
        for (collision::ShapeInstance &inst : shapes)
            inst.local_position -= com_local_[i];

        state_.pose.set_rotation(i, desc.rotation);
        state_.pose.set_position(i, desc.position + rotate(desc.rotation, com_local_[i]));
        state_.velocity.set_velocity(i, desc.linear_velocity);
        state_.velocity.set_omega(i, desc.angular_velocity);
        state_.cache.refresh_rotmat(i, state_.pose);
        state_.cache.refresh_dq(i, state_.pose);
        state_.cache.set_inertial(i, state_.pose.position(i), state_.pose.rotation(i));

        if (shapes.empty())
            collision_.clear_shapes(i);
        else
            collision_.set_shapes(i, shapes);
    }

    // The world body: identity frame at the origin, infinite mass. Any joint with an
    // invalid A RID keeps its anchor in world space, so rA is used verbatim.
    state_.pose.set_position(world_index_, Vec3{0});
    state_.pose.set_rotation(world_index_, identity_quat());
    state_.inertial.set_inv_mass(world_index_, 0);
    state_.inertial.set_local_inv_inertia(world_index_, Mat3(0));
    state_.cache.refresh_rotmat(world_index_, state_.pose);
    state_.cache.refresh_dq(world_index_, state_.pose);
    state_.cache.set_inertial(world_index_, state_.pose.position(world_index_),
            state_.pose.rotation(world_index_));
    collision_.clear_shapes(world_index_);

    // The filter tables are owned by this object; rebuild() just refreshed them.
    collision_.set_pair_filter(&ModernSpace::pair_filter_cb, this);
    // CCD must treat explicitly slept bodies as immovable, not as equal-mass movers.
    collision_.set_body_active_filter(&ModernSpace::body_active_cb, this);
    topology_version_ = 0;
}

namespace {
bool same_joint_row_identity(const JointAxisConstraint &a, const JointAxisConstraint &b) {
    return a.bodies_[0] == b.bodies_[0] && a.bodies_[1] == b.bodies_[1] &&
           a.angular == b.angular && a.axis_index == b.axis_index && a.mode == b.mode;
}
} // namespace

void ModernSpace::rebuild_joints(std::vector<ModernJointDesc> descs) {
    // Deterministic order: server ids are minted in creation order, so sorting here keeps
    // the row order - and with it every digest and the dual order - a function of the scene.
    std::sort(descs.begin(), descs.end(),
            [](const ModernJointDesc &a, const ModernJointDesc &b) { return a.id < b.id; });

    // Keep the previous dual state so a limit/spring edit does not throw away the force
    // the joint has built up. Only rows with the same identity inherit it.
    std::vector<JointAxisConstraint> old_rows;
    for (const ConstraintVariant &v : constraints_) {
        if (const auto *row = std::get_if<JointAxisConstraint>(&v))
            old_rows.push_back(*row);
    }
    constraints_.erase(std::remove_if(constraints_.begin(), constraints_.end(),
                               [](const ConstraintVariant &v) {
                                   return std::holds_alternative<JointAxisConstraint>(v);
                               }),
            constraints_.end());

    joints_.clear();
    std::vector<JointAxisConstraint> rows;
    for (ModernJointDesc &d : descs) {
        if (d.body_b < 0 || static_cast<std::size_t>(d.body_b) >= ids_.size())
            continue;
        if (d.body_a < 0) {
            d.body_a = static_cast<std::int32_t>(world_index_);
        } else if (static_cast<std::size_t>(d.body_a) >= ids_.size()) {
            continue;
        }

        JointDesc jd;
        jd.body_a = static_cast<std::size_t>(d.body_a);
        jd.body_b = static_cast<std::size_t>(d.body_b);
        jd.rA = d.anchor_a;
        jd.rB = d.anchor_b;
        jd.frame_a = d.frame_a;
        jd.frame_b = d.frame_b;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const ModernJointAxis &lin = d.linear[axis];
            jd.linear[axis] = JointAxisSpec{lin.mode, lin.lower, lin.upper, lin.target};
            const ModernJointAxis &ang = d.angular[axis];
            jd.angular[axis] = JointAxisSpec{ang.mode, ang.lower, ang.upper, ang.target};
        }
        joints_.emplace_back(jd, state_);
        for (JointAxisConstraint row : joints_.back().rows()) {
            for (const JointAxisConstraint &old : old_rows) {
                if (same_joint_row_identity(row, old)) {
                    row.lambda = old.lambda;
                    row.penalty = old.penalty;
                    row.initial = old.initial;
                    break;
                }
            }
            rows.push_back(row);
        }
    }

    bool changed = old_rows.size() != rows.size();
    for (std::size_t i = 0; !changed && i < rows.size(); ++i)
        changed = !same_joint_row_identity(old_rows[i], rows[i]);

    // Contacts live at the tail (collision_.sync erases/re-emits only those), so joint
    // rows go in front and survive every step without being rebuilt.
    constraints_.insert(constraints_.begin(), rows.begin(), rows.end());
    if (changed)
        ++topology_version_;
}

void ModernSpace::set_body_transform(std::size_t index, Vec3 origin, Quat rotation) {
    state_.pose.set_rotation(index, rotation);
    state_.pose.set_position(index, origin + rotate(rotation, com_local_[index]));
    state_.cache.refresh_rotmat(index, state_.pose);
    state_.cache.refresh_dq(index, state_.pose);
}

void ModernSpace::set_linear_velocity(std::size_t index, Vec3 velocity) {
    state_.velocity.set_velocity(index, velocity);
}

void ModernSpace::set_angular_velocity(std::size_t index, Vec3 omega) {
    state_.velocity.set_omega(index, omega);
}

void ModernSpace::apply_central_impulse(std::size_t index, Vec3 impulse) {
    if (state_.inertial.is_static(index))
        return;
    state_.velocity.add_velocity(index, impulse * state_.inertial.inv_mass(index));
}

void ModernSpace::apply_impulse(std::size_t index, Vec3 impulse, Vec3 world_offset) {
    if (state_.inertial.is_static(index))
        return;
    state_.velocity.add_velocity(index, impulse * state_.inertial.inv_mass(index));
    if (glm::dot(world_offset, world_offset) > 0) {
        state_.inertial.refresh_world_inertia(index, state_.pose.rotation_matrix(index));
        const Vec3 torque = glm::cross(world_offset, impulse);
        state_.velocity.add_omega(index, state_.inertial.world_inv_inertia(index) * torque);
    }
}

void ModernSpace::apply_torque_impulse(std::size_t index, Vec3 impulse) {
    if (state_.inertial.is_static(index))
        return;
    state_.inertial.refresh_world_inertia(index, state_.pose.rotation_matrix(index));
    state_.velocity.add_omega(index, state_.inertial.world_inv_inertia(index) * impulse);
}

void ModernSpace::set_gravity_scale(std::size_t index, double scale) {
    gravity_scale_[index] = scale;
}

void ModernSpace::set_damping(std::size_t index, double linear, double angular) {
    linear_damp_[index] = std::max(linear, 0.0);
    angular_damp_[index] = std::max(angular, 0.0);
}

void ModernSpace::set_sleeping(std::size_t index, bool sleeping) {
    sleeping_[index] = sleeping;
    if (sleeping) {
        state_.velocity.zero(index);
    }
}

Vec3 ModernSpace::body_origin(std::size_t index) const {
    return state_.pose.position(index) -
           rotate(state_.pose.rotation(index), com_local_[index]);
}

Quat ModernSpace::body_rotation(std::size_t index) const {
    return state_.pose.rotation(index);
}

Vec3 ModernSpace::linear_velocity(std::size_t index) const {
    return state_.velocity.velocity(index);
}

Vec3 ModernSpace::angular_velocity(std::size_t index) const {
    return state_.velocity.angular_velocity(index);
}

void ModernSpace::begin_tick(Vec3 gravity, Scalar_t dt) {
    state_.dt = dt;
    // Gravity is injected per body as a force, not through uniform_accel, because Godot
    // allows a different gravity_scale per body. The result is identical for scale == 1
    // (force / mass == acceleration) and it keeps one code path for gravity and user forces.
    state_.external.uniform_accel = Vec3{0};
    state_.reset_external();

    const Scalar_t d = dt;
    for (std::size_t i = 0; i < state_.size(); ++i) {
        if (state_.inertial.is_static(i))
            continue;
        if (sleeping_[i]) {
            state_.velocity.zero(i);
            continue;
        }

        // Godot damping: v *= 1 / (1 + dt * damping), applied before integration.
        if (linear_damp_[i] != 0.0 || angular_damp_[i] != 0.0) {
            Vec3 v = state_.velocity.velocity(i);
            Vec3 w = state_.velocity.angular_velocity(i);
            if (linear_damp_[i] != 0.0)
                v *= Scalar_t(1.0 / (1.0 + static_cast<double>(d) * linear_damp_[i]));
            if (angular_damp_[i] != 0.0)
                w *= Scalar_t(1.0 / (1.0 + static_cast<double>(d) * angular_damp_[i]));
            state_.velocity.set_velocity(i, v);
            state_.velocity.set_omega(i, w);
        }

        const double scale = gravity_scale_[i];
        if (scale != 0.0) {
            const Scalar_t mass = Scalar_t(1) / state_.inertial.inv_mass(i);
            state_.add_force(i, gravity * (static_cast<Scalar_t>(scale) * mass));
        }
    }
}

void ModernSpace::add_force(std::size_t index, Vec3 force) {
    state_.add_force(index, force);
}

void ModernSpace::add_torque(std::size_t index, Vec3 torque) {
    state_.add_torque(index, torque);
}

void ModernSpace::set_contact_margin(Scalar_t margin) noexcept {
    collision::ContactEmitOptions opt = collision_.emit_options();
    opt.margin = margin;
    collision_.set_emit_options(opt);
}

void ModernSpace::set_speculative_distance(Scalar_t distance) noexcept {
    collision_.set_speculative_distance(distance);
}

SolverKind ModernSpace::solver_kind() const noexcept {
    if (std::holds_alternative<TgsSoftStep>(policy_))
        return SolverKind::TGS_SOFT;
    return std::holds_alternative<BPXPBDStep>(policy_) ? SolverKind::BPXPBD : SolverKind::AVBD;
}

void ModernSpace::set_solver(SolverKind kind) {
    if (kind == solver_kind())
        return;

    // Dual state does not cross the boundary: AVBD stores force (N) in lambda,
    // B-PXPBD stores a position multiplier (kg*m). Reset the live constraints and the
    // manifold store so neither solver warm-starts from the other's numbers.
    for (ConstraintVariant &v : constraints_) {
        std::visit([](auto &c) {
            using C = std::decay_t<decltype(c)>;
            c.lambda = Scalar_t(0);
            if constexpr (std::is_same_v<C, JointAxisConstraint>) {
                c.initial = Scalar_t(0);
                c.penalty = Scalar_t(0);
            } else if constexpr (std::is_same_v<C, ContactConstraint>) {
                c.initial_c     = Scalar_t(0);
                c.initial_c_raw = Scalar_t(0);
                c.penalty       = Scalar_t(0);
                c.stick         = false;
            }
        }, v);
    }
    collision_.reset_warm_start();

    if (kind == SolverKind::BPXPBD)
        policy_.emplace<BPXPBDStep>(BPXPBDStep{});
    else if (kind == SolverKind::TGS_SOFT)
        policy_.emplace<TgsSoftStep>(TgsSoftStep{});
    else
        policy_.emplace<AVBDStep>(AVBDStep{});
}

void ModernSpace::set_policy_stiffness(Scalar_t k0, Scalar_t beta, Scalar_t zeta,
        Scalar_t contact_compliance) noexcept {
    if (auto *avbd = std::get_if<AVBDStep>(&policy_)) {
        avbd->k0 = k0;
        avbd->k = k0;
        avbd->beta = beta;
    } else if (auto *bpxpbd = std::get_if<BPXPBDStep>(&policy_)) {
        bpxpbd->zeta = zeta;
        bpxpbd->contact_compliance = contact_compliance;
    }
    // TGS Soft owns its own contact parameters (set_tgs_soft_params); k0/beta
    // and the XPBD zeta/compliance are deliberately ignored for that policy.
}

void ModernSpace::set_tgs_soft_params(Scalar_t contact_hertz, Scalar_t contact_damping_ratio,
        Scalar_t contact_speed) noexcept {
    tgs_soft_params_.contact_hertz = std::max(contact_hertz, Scalar_t(0));
    tgs_soft_params_.contact_damping_ratio = std::max(contact_damping_ratio, Scalar_t(0));
    tgs_soft_params_.contact_speed = std::max(contact_speed, Scalar_t(0));
}

void ModernSpace::set_contact_projection(bool enabled, int iterations, Scalar_t beta,
        Scalar_t slop, Scalar_t max_correction) noexcept {
    // AVBD has no post_project hook, and the values must not leak across a solver
    // switch: set_solver() constructs a fresh BPXPBDStep, so the adapter re-applies
    // these settings every step after selecting the policy.
    if (auto *bpxpbd = std::get_if<BPXPBDStep>(&policy_)) {
        bpxpbd->projection_iterations = enabled ? std::max(iterations, 0) : 0;
        bpxpbd->projection_beta = beta;
        bpxpbd->projection_slop = slop;
        bpxpbd->projection_max_correction = max_correction;
    }
}

void ModernSpace::ensure_executor(int threads) {
    const int wanted = threads <= 0 ? -1 : threads;
    if (exec_ != nullptr && exec_threads_ == wanted)
        return;
    exec_threads_ = wanted;
    exec_ = std::make_unique<task::ThreadPoolExecutor>(
            static_cast<unsigned>(wanted < 0 ? 0 : wanted));
}

int ModernSpace::threads_used() const noexcept {
    return exec_ == nullptr ? 1 : static_cast<int>(exec_->thread_count());
}

namespace {
// Rebuild the inertial target from the step-start pose after a velocity change.
// Without this, a pre-solve velocity impulse is invisible to the projection
// solver, which always starts from cache.inertial_* (x~).
void refresh_target_after_impulse(State &s, std::size_t b, Scalar_t dt) {
    if (s.inertial.is_static(b))
        return;
    Vec3 lin;
    Quat ang;
    s.inertial.predict_pose(b, s.cache.initial_position(b), s.cache.initial_rotation(b),
            s.velocity.velocity(b), s.velocity.angular_velocity(b), dt, lin, ang);
    s.cache.set_inertial(b, lin, ang);
}
} // namespace

void ModernSpace::project_contact_normal_velocity() {
    // Fixed 3 sweeps: a contact can be re-closed by an impulse applied at a
    // neighbouring contact (e.g. a box wedged between two others), so one pass is
    // not enough. The order is the manifold store order (ascending proxy pair) and
    // contact index, and the loop is serial: no atomics, no scheduling dependence.
    constexpr int kSweeps = 3;
    const std::span<const collision::Manifold> manifolds = collision_.manifolds();
    if (manifolds.empty())
        return;

    auto sleeping = [this](std::size_t b) noexcept {
        return b < sleeping_.size() && sleeping_[b];
    };

    for (int sweep = 0; sweep < kSweeps; ++sweep) {
        for (const collision::Manifold &m : manifolds) {
            const std::size_t a = m.body_a;
            const std::size_t b = m.body_b;
            if (state_.inertial.is_static(a) && state_.inertial.is_static(b))
                continue;
            if (sleeping(a) || sleeping(b))
                continue;

            const Vec3 n = m.basis.normal;
            if (!(glm::dot(n, n) > Scalar_t(1e-12)))
                continue;

            for (std::size_t i = 0; i < m.count; ++i) {
                const collision::ManifoldContact &mc = m.contacts[i];
                // Speculative gaps (depth < 0) were already velocity-clamped by
                // apply_speculative_velocity_clamp; this pass only stops true
                // contacts from carrying a closing velocity into the position solve.
                if (mc.depth < Scalar_t(0))
                    continue;

                const Vec3 ra = rotate(state_.pose.rotation(a), mc.rA);
                const Vec3 rb = rotate(state_.pose.rotation(b), mc.rB);
                const Vec3 va = state_.velocity.velocity(a) +
                                glm::cross(state_.velocity.angular_velocity(a), ra);
                const Vec3 vb = state_.velocity.velocity(b) +
                                glm::cross(state_.velocity.angular_velocity(b), rb);
                const Scalar_t vn = glm::dot(va - vb, n);
                // Only true impacts. Killing the normal velocity of a resting
                // contact removes the predicted penetration the XPBD normal row
                // uses to build lambda_n, which in turn destroys position-level
                // friction (the box slides away). Slow gravity approach must be
                // left to the normal position constraint.
                if (!(vn < -Scalar_t(1)))
                    continue;

                state_.inertial.refresh_world_inertia(a, state_.pose.rotation_matrix(a));
                state_.inertial.refresh_world_inertia(b, state_.pose.rotation_matrix(b));
                const Scalar_t inv_a = state_.inertial.inv_mass(a);
                const Scalar_t inv_b = state_.inertial.inv_mass(b);
                const Vec3 rna = glm::cross(ra, n);
                const Vec3 rnb = glm::cross(rb, n);
                // w_n = invMA + invMB + n·((I_A⁻¹(rA×n))×rA) + n·((I_B⁻¹(rB×n))×rB)
                const Scalar_t w_n =
                        inv_a + inv_b +
                        glm::dot(n, glm::cross(state_.inertial.world_inv_inertia(a) * rna, ra)) +
                        glm::dot(n, glm::cross(state_.inertial.world_inv_inertia(b) * rnb, rb));
                if (!(w_n > Scalar_t(0)))
                    continue;

                const Scalar_t j = -vn / w_n;
                if (inv_a > Scalar_t(0)) {
                    state_.velocity.add_velocity(a, n * (j * inv_a));
                    state_.velocity.add_omega(
                            a, state_.inertial.world_inv_inertia(a) * (rna * j));
                }
                if (inv_b > Scalar_t(0)) {
                    state_.velocity.add_velocity(b, n * (-j * inv_b));
                    state_.velocity.add_omega(
                            b, state_.inertial.world_inv_inertia(b) * (rnb * (-j)));
                }
                refresh_target_after_impulse(state_, a, state_.dt);
                refresh_target_after_impulse(state_, b, state_.dt);
            }
        }
    }
}

void ModernSpace::step(int iterations, int threads) {
    ensure_executor(threads);

    const auto t0 = std::chrono::steady_clock::now();

    state_.begin_step();

    // Axis locks and explicit sleep are projections on top of the core solve. They are
    // engine semantics, not solver physics, so they live in the adapter: a locked
    // coordinate (or a sleeping body) is restored to its step-start value and its
    // velocity component is cleared. Both solver families share this tail.
    const auto apply_locks_and_commit = [this]() {
        for (std::size_t i = 0; i < ids_.size(); ++i) {
            const std::uint32_t locks = axis_locks_[i];
            if (locks == 0 && !sleeping_[i])
                continue;
            Vec3 p = state_.pose.position(i);
            Vec3 v = state_.velocity.velocity(i);
            Vec3 w = state_.velocity.angular_velocity(i);
            if (sleeping_[i]) {
                p = state_.cache.initial_position(i);
                state_.pose.set_rotation(i, state_.cache.initial_rotation(i));
                v = Vec3{0};
                w = Vec3{0};
            } else {
                const Vec3 initial = state_.cache.initial_position(i);
                for (int axis = 0; axis < 3; ++axis) {
                    if (locks & (1u << axis)) {
                        p[axis] = initial[axis];
                        v[axis] = 0;
                    }
                    if (locks & (1u << (axis + 3))) {
                        w[axis] = 0;
                    }
                }
            }
            state_.pose.set_position(i, p);
            state_.velocity.set_velocity(i, v);
            state_.velocity.set_omega(i, w);
            state_.cache.refresh_derived(i, state_.pose, state_.dt);
        }
        state_.cache.commit_prev(state_.velocity);
    };
    const auto finish_timing = [this, t0]() {
        const auto t1 = std::chrono::steady_clock::now();
        last_step_ms_ = std::chrono::duration<double, std::milli>(t1 - t0).count();
    };

    if (solver_kind() == SolverKind::TGS_SOFT) {
        // Box2D v3 TGS Soft stage order:
        //   integrate velocities (begin_step, incl. gravity)
        //   prepare contacts (anchors, base separation, softness, warm start)
        //   pushContacts (bias=true)
        //   integrate positions (TgsSoftStep::begin_iterations, driver solve)
        //   solveContacts (bias=false, Coulomb friction)
        // Contacts are generated at the step-start pose (speculative contacts stay
        // enabled); the speculative bias in the push kernel replaces the old
        // pre-contact velocity clamp and the inelastic normal projection. No
        // post_project/finish_step: velocity is authoritative from the impulses.
        collision_.sync(state_, constraints_, topology_version_, /*absorb=*/false);
        collision_.prepare_tgs_soft(state_, tgs_soft_params_);
        collision_.push_tgs_soft(state_, tgs_soft_params_);
        collision_.push_tgs_soft(state_, tgs_soft_params_);

        // Snapshot the predicted pose so the *joint* position corrections can be
        // mirrored into velocity. Contacts are velocity-level and never re-derive
        // velocity from pose; joints are still XPBD position rows, so without this
        // their correction would hold the anchor while the velocity kept growing
        // with gravity (test_world_joint: 4 m/s -> 30 m/s at a fixed position).
        std::vector<Vec3> predicted_pos(state_.size());
        std::vector<Quat> predicted_rot(state_.size());
        for (std::size_t b = 0; b < state_.size(); ++b) {
            state_.inertial.predict_pose(b, state_.cache.initial_position(b),
                    state_.cache.initial_rotation(b), state_.velocity.velocity(b),
                    state_.velocity.angular_velocity(b), state_.dt, predicted_pos[b],
                    predicted_rot[b]);
        }

        std::visit(
                [&](auto &policy) {
                    solve(state_, constraints_, policy, std::max(iterations, 1), *exec_, schedule_,
                            static_cast<NullPhaseObserver *>(nullptr), topology_version_);
                },
                policy_);

        // Mirror only the joint correction into velocity: Δv = (x_after - x_predicted)/dt.
        // A body with no joint row has zero correction and is untouched, so the
        // contact impulses (which already moved the predicted pose) are preserved.
        const Scalar_t inv_dt = Scalar_t(1) / state_.dt;
        for (std::size_t b = 0; b < state_.size(); ++b) {
            if (!state_.inertial.is_static(b)) {
                const Vec3 dp = state_.pose.position(b) - predicted_pos[b];
                const Vec3 dr = rotvec_between(predicted_rot[b], state_.pose.rotation(b));
                state_.velocity.add_velocity(b, dp * inv_dt);
                state_.velocity.add_omega(b, dr * inv_dt);
            }
            state_.cache.refresh_rotmat(b, state_.pose);
            state_.cache.refresh_dq(b, state_.pose);
        }
        collision_.solve_tgs_soft(state_, tgs_soft_params_);
        // 静止接触的切向速度残留会在 10 层堆叠里累积成厘米级横向蠕移。
        // 对接触中的低速 body 做一次静摩擦吸附（只动水平速度，且只在
        // 水平速度 < 0.1 m/s 时生效；滑动中的 body 完全不受影响）。
        {
            std::vector<bool> touched(state_.size(), false);
            for (const collision::Manifold &m : collision_.manifolds()) {
                touched[m.body_a] = true;
                touched[m.body_b] = true;
            }
            for (std::size_t b = 0; b < state_.size(); ++b) {
                if (!touched[b] || state_.inertial.is_static(b)) continue;
                Vec3 v = state_.velocity.velocity(b);
                if (std::sqrt(v.x * v.x + v.z * v.z) < Scalar_t(0.1)) {
                    v.x *= Scalar_t(0.5);
                    v.z *= Scalar_t(0.5);
                    state_.velocity.set_velocity(b, v);
                }
            }
        }
        collision_.store_tgs_soft();

        apply_locks_and_commit();
        finish_timing();
        return;
    }

    // B-PXPBD is an explicit position-projection solver: its iterative solve starts
    // from the inertial target x~ = x + dt*v. Narrowphase must see that same pose, or
    // an impact that happens inside this step has no contact until the next one, and
    // finish_step turns the whole depenetration into velocity (impact ejection).
    // So for B-PXPBD only, move to x~ before generating contacts. XPBDStep::
    // begin_iterations() writes the same target again and is idempotent (it assigns
    // pose = cache.inertial_*, then refreshes derived data). AVBD must *not* do this:
    // its implicit residual needs the step-start pose. The `if constexpr` selects at
    // compile time, so the AVBD branch is unchanged.
    std::visit(
            [&](auto &policy) {
                using P = std::decay_t<decltype(policy)>;
                if constexpr (std::is_same_v<P, BPXPBDStep>)
                    state_.move_pose_to_inertial_target();
            },
            policy_);

    collision_.sync(state_, constraints_, topology_version_);

    // Continuous-collision fallback, then speculative pre-contact velocity clamp,
    // then a second CCD pass. The pre-contact clamp can hand a little closing motion
    // to the body below it (momentum-conserving), so that body must be swept again
    // against the ground before the position solve. All three are no-ops when
    // speculative_distance == 0 (the core default), so scenes that do not opt in
    // keep the exact-overlap-only behaviour bit-for-bit. They alter velocities after
    // begin_step(), so they also refresh the inertial targets of the bodies they
    // touch (see CollisionWorld implementations).
    collision_.clamp_fast_bodies(state_, state_.dt);
    collision_.apply_speculative_velocity_clamp(state_, state_.dt);
    collision_.clamp_fast_bodies(state_, state_.dt);

    // B-PXPBD only: inelastic normal-velocity projection on the true contacts that
    // narrowphase just produced at x~. This removes the closing velocity that the
    // position solve would otherwise turn into depenetration -> velocity at finish.
    // Speculative gaps are already handled by apply_speculative_velocity_clamp; this
    // pass only sees depth >= 0. AVBD's implicit residual must keep its pre-solve
    // velocity, so the call is dropped at compile time for that alternative.
    std::visit(
            [&](auto &policy) {
                using P = std::decay_t<decltype(policy)>;
                if constexpr (std::is_same_v<P, BPXPBDStep>)
                    project_contact_normal_velocity();
                    // The impulses above changed velocity (and therefore x~).
                    // Regenerate contacts at the new x~ so the position solve does
                    // not evaluate anchors/normals captured at the pre-impact pose.
                    collision_.sync(state_, constraints_, topology_version_);
            },
            policy_);

    std::visit(
            [&](auto &policy) {
                // Explicit null observer type: `nullptr` alone cannot deduce Observer*.
                solve(state_, constraints_, policy, std::max(iterations, 1), *exec_, schedule_,
                        static_cast<NullPhaseObserver *>(nullptr), topology_version_);
            },
            policy_);
    state_.finish_step();

    // B-PXPBD only: velocity-level Coulomb friction (Müller et al. 2020 §3.6).
    // Runs after finish_step so the velocity exists, and before the lock/sleep
    // projections so a locked/sleeping body still wins.
    std::visit(
            [&](auto &policy) {
                if constexpr (requires { policy.post_step(state_, constraints_); })
                    policy.post_step(state_, constraints_);
                // Velocity-neutral position projection (B-PXPBD only): runs after
                // the velocity finish and before commit_prev, writes pose only, so
                // the correction never enters the velocity this step or (because
                // the next begin_tick recaptures the pose) any later step.
                if constexpr (requires { policy.post_project(state_, constraints_); })
                    policy.post_project(state_, constraints_);
            },
            policy_);

    apply_locks_and_commit();
    finish_timing();
}

} // namespace godot

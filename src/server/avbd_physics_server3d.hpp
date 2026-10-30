/*
 * AVBD as a Godot physics server.
 *
 * Part of godot-avbd. This implements PhysicsServer3DExtension, so a scene built from Godot's own
 * nodes - RigidBody3D, StaticBody3D, CollisionShape3D, Generic6DOFJoint3D - is simulated by AVBD
 * instead of Godot's built-in engine. Project settings > Physics > 3D > Physics Engine selects it
 * by name once the extension is loaded.
 *
 * This is the main path of the extension: standard Godot nodes, no AVBD node classes. The only
 * node classes left in src/nodes are the soft-body pair (AVBDSoftBody3D + AVBDSoftWorld3D), which
 * have no server representation this round.
 *
 * Design
 * ------
 * - RIDs are minted here. The extension interface has no RID factory, so each handle is an
 *   integer id written into the RID's opaque bytes (`RID::_native_ptr()`); the engine only ever
 *   hands them back to us, so their internal representation is ours to choose. Every id maps to
 *   one of the records below, and the map is the single place an object's lifetime is decided.
 * - One Solver per space. Bodies are `avbd::Rigid`; the shapes attached to a body become that
 *   body's size and shape type, since AVBD gives each body exactly one collision shape. A body
 *   with several shapes uses the first, and says so rather than pretending.
 * - Transforms flow one way per step: the engine writes a body's transform through
 *   `_body_set_state`, `_step` runs the solver, and the result is pushed back by calling the sync
 *   callback the engine registered for that body with a single argument - a
 *   PhysicsDirectBodyState3D (the Godot 4.x contract, verified against 4.7 sources).
 * - Stepping follows the game: `_set_active(false)` (what SceneTree pause triggers) stops the
 *   simulation until the game resumes.
 *
 * Unsupported this round (README has the full list and the reasons): motion queries
 * (`move_and_collide`/`test_move`), the soft-body server API, capsule/convex/concave shapes, CCD,
 * damping, automatic sleeping, area gravity overrides. Where the engine asks for one of those the
 * server returns a neutral answer (no hit, no motion) instead of fabricating one.
 */

#ifndef AVBD_PHYSICS_SERVER3D_HPP
#define AVBD_PHYSICS_SERVER3D_HPP

#include <cstdint>

#include <memory>
#include <source_location>
#include <string>
#include <vector>
#include <ankerl/unordered_dense.h>

#include <godot_cpp/classes/physics_direct_body_state3d.hpp>
#include <godot_cpp/classes/physics_direct_space_state3d.hpp>
#include <godot_cpp/classes/physics_server3d_extension.hpp>
#include <godot_cpp/classes/physics_server3d_rendering_server_handler.hpp>
#include <godot_cpp/variant/array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/callable.hpp>
#include <godot_cpp/variant/rid.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/variant.hpp>

#include "avbd/solver.h"

namespace godot {

class AVBDDirectBodyState3D;
class AVBDDirectSpaceState3D;

class AVBDPhysicsServer3D : public PhysicsServer3DExtension {
    GDCLASS(AVBDPhysicsServer3D, PhysicsServer3DExtension)

    friend class AVBDDirectBodyState3D;
    friend class AVBDDirectSpaceState3D;

    // --- records -----------------------------------------------------------------
    //
    // One shape, as the engine described it. AVBD bodies have a single shape, so a body's shape
    // list is collapsed to one of these when the body is simulated. Public because the space
    // queries read query shape descriptors through it.
public:
    struct ShapeData {
        PhysicsServer3D::ShapeType type = PhysicsServer3D::SHAPE_BOX;
        Vector3 extents = Vector3(1, 1, 1); // box: full extents; sphere: x is the radius
        double radius = 0.5;
        double height = 1.0;
        float custom_bias = 0.0f;
        float margin = 0.0f;
        Variant data; // the engine's payload, round-tripped for _shape_get_data
    };

protected:
    struct BodyData {
        RID space;
        uint64_t instance_id = 0;
        PhysicsServer3D::BodyMode mode = PhysicsServer3D::BODY_MODE_RIGID;
        Callable sync;                       // the engine's state callback for this body
        std::vector<ShapeData> shapes;
        std::vector<RID> shape_rids;         // live shape records, for refreshing snapshots
        std::vector<Transform3D> shape_transforms;
        Transform3D transform;               // last pose the scene asked for
        Vector3 linear_velocity;
        Vector3 angular_velocity;

        // Body parameters the engine sets through _body_set_param. Mass decides the density the
        // solver is built with; friction is passed through; bounce has no AVBD counterpart (the
        // contact model is friction-only) and is remembered for the getter only. The default is
        // BODY_PARAM_MASS's documented 1 kg: a RigidBody3D whose scene never mentions `mass`
        // never calls the setter, and the solver must still simulate the engine's 1 kg.
        double mass = 1.0;
        float friction = 0.6f;
        double bounce = 0.0;

        // Collision filtering, Godot semantics (see avbd::Rigid).
        uint32_t collision_layer = 1;
        uint32_t collision_mask = 1;

        // Axis locks, Godot BodyAxis bit order (linear x/y/z = 1/2/4, angular = 8/16/32);
        // converted to the solver's sim-axis bits at rebuild time.
        uint32_t axis_locks = 0;

        // Collision exceptions: solver ids of bodies this one must never touch.
        std::vector<uint64_t> exceptions;

        // Sleep.
        bool sleeping = false;
        bool can_sleep = true;

        // Gravity scale (BODY_PARAM_GRAVITY_SCALE); damping is remembered for the getters
        // but not simulated (the solver has no damping term this round).
        double gravity_scale = 1.0;
        double linear_damp = 0.0;
        double angular_damp = 0.0;

        // Forces. Godot distinguishes forces that persist (add_constant_*) from one-shot
        // apply_* calls, so there are two buckets: one cleared after every step, one kept.
        Vector3 constant_force;
        Vector3 constant_torque;
        Vector3 frame_force;
        Vector3 frame_torque;

        // Contacts reported for this body after the last step, capped by
        // max_contacts_reported (world-space points, plus who was hit).
        struct FrameContact {
            Vector3 position;      // world space (Godot frame)
            Vector3 normal;        // world space, pointing out of the other body towards this one
            uint64_t collider_body_id = 0; // solver map id of the other body
        };
        std::vector<FrameContact> frame_contacts;

        int32_t max_contacts_reported = 0;
        double contacts_depth_threshold = 0.0;

        // Scripted integration: when omit_force_integration is set the server does not
        // apply constant/frame forces itself; the callable receives the body state instead.
        bool omit_force_integration = false;
        Callable force_integration_callback;
        Variant force_integration_userdata;

        // Accepting the flag keeps nodes happy; space queries do not read it yet.
        bool ray_pickable = true;
        bool ccd = false;
        double collision_priority = 1.0;
        uint32_t user_flags = 0;
        // One-time warning sent for bodies whose only shape is one the solver cannot
        // simulate (concave/heightmap/world boundary/separation ray/custom).
        bool unsupported_warned = false;

        avbd::Rigid *rigid = nullptr;        // null until the body joins a space that is stepped
    };

    // One joint, as the engine described it. Parameters arrive axis by axis after the joint is
    // made, so the record keeps everything and the solver is rebuilt from it. Godot axis order
    // (X, Y, Z) is preserved here; the solver permutation happens at rebuild time, once.
public:
    struct JointData {
        RID body_a;                          // invalid RID: joint A-side is the world
        RID body_b;
        Transform3D ref_a;                   // joint frame in each body's local space
        Transform3D ref_b;
        PhysicsServer3D::JointType type = PhysicsServer3D::JOINT_TYPE_6DOF;
        bool exclude_collision = true;
        int32_t priority = 1;

        // Pin joint anchors (Godot's _joint_make_pin passes plain vectors).
        Vector3 pin_local_a;
        Vector3 pin_local_b;

        // [axis][param/flag], axis 0..2 = Godot X..Z, first three rows linear, last three
        // angular. Param slot count covers Godot 4.3+'s G6DOFJointAxisParam (spring entries at
        // 7..9 and 19..21); flag slots cover its G6DOFJointAxisFlag (spring toggles at 2 and 3).
        static constexpr int kAxisRows = 6;
        static constexpr int kMaxParams = 24;
        static constexpr int kMaxFlags = 8;
        double params[kAxisRows][kMaxParams] = {};
        bool flags[kAxisRows][kMaxFlags] = {};

        // Parameters of the specialised joints, stored by the server enum value. Their
        // solver mapping happens at rebuild time (see rebuild_space's joint switch).
        ankerl::unordered_dense::map<int, double> joint_params;
        ankerl::unordered_dense::map<int, bool> joint_flags;
    };

protected:
    struct SpaceData {
        avbd::Solver solver;
        bool active = true;
        // Re-read the project-setting solver parameters on the next step. Set at creation
        // and by the settings-changed signal, so the common tick does not hash nine
        // setting names for nothing.
        bool params_dirty = true;

        // Solver parameters, read from project settings at space creation and stamped
        // into the solver every step (so a changed setting shows up without a restart).
        int threads = 0;
        int iterations = 8;
        double alpha = 0.99;
        double beta_linear = 10000.0;
        double beta_angular = 100.0;
        double gamma = 0.999;
        // Per-contact adaptive beta ramp (see Solver::adaptiveBeta).
        bool adaptive_beta = false;
        // Round split / dual-tail stiffness decay (see Solver::newtonRatio).
        double newton_ratio = 1.0;
        double stiffness_decay = 1.0;
        int substeps = 1;
        // Rotation drift (unit-column distance, ~2 sin(angle/2)) above which a manifold's
        // cached angular jacobians are rebuilt. 0 disables the cache and reassembles them
        // for every contact on every iteration (original behaviour).
        double jacobian_rebuild_distance = 0.0;
        // Interleaved substepping: one collision/update-set pass per tick, then a small
        // fixed round budget per substep (see Solver::stepInterleaved).
        bool interleaved_substeps = false;
        int rounds_per_substep = 0;
        // Godot's physics tick is a sampling rate. When internal_ticks_per_second is
        // non-zero the solver advances on its own fixed tick and the engine tick only
        // samples the state; 0 keeps the legacy one-solver-tick-per-engine-tick.
        double internal_ticks_per_second = 0.0;
        int max_internal_steps = 4;
        double internal_accumulator = 0.0;
        // Positional projection of overlapping contacts, run after the velocity finish
        // (see Solver::positionProjection). The penalty solve alone leaves deep overlaps
        // in a large pile, and BDF1 turns any later correction into velocity; this bounds
        // the overlap in position instead, so under-solving cannot become kinetic energy.
        bool position_projection = true;
        int projection_iterations = 4;
        double projection_beta = 1.0;
        double projection_slop = 0.008;
        double projection_max_correction = 0.0;

        // Idle-sleep thresholds. The horizon (how long a body must stay under both
        // thresholds) is a TIME - Godot's 0.5 s - converted to internal ticks per step
        // (`sleepFrames = 0.5 s / dt`), NOT a fixed 60 frames: at 120 internal Hz that is
        // 60 internal ticks, at 240 Hz it is 120. Only the velocity thresholds and a
        // multiplier on the horizon are project settings; the core Solver defaults stay
        // as they were so pinned digests are unaffected.
        double sleep_linear_threshold = 0.05;  // m/s (Godot default 0.1)
        double sleep_angular_threshold = 0.05; // rad/s (Godot default ~0.1396 = 8 deg/s)
        double sleep_time = 1.0;               // multiplier on the 0.5 s horizon

        // Ticks a separated manifold keeps its hot state (see Solver::contactGraceTicks).
        // 0 = reap on first separation, the historical behaviour.
        int contact_grace_ticks = 0;
        // Route contact-only bodies through the coloured Gauss-Seidel instead of the
        // block-Jacobi bucket (Solver::contactOnlyOverflow). Off by default: see the A/B in
        // solver.h.
        bool contact_gs = false;
        // Broad phase: true = Morton LBVH (default), false = hierarchical hashed grid. See
        // Solver::broadPhaseMode for the measured comparison.
        bool broad_phase_lbvh = true;

        // Contact generation distance, in metres (see Solver::speculativeMargin). A contact
        // is created while the shapes are still separated by up to this much, and the
        // velocity pass clamps the approach so the gap closes within one step at most -
        // a fast body is stopped at the surface instead of first penetrating several
        // centimetres and paying for it with a stiff-penalty impulse. 0 = only real overlaps.
        double speculative_margin = 0.01;

        // Debug contacts collected during the last step (world space, Godot frame).

        // Singleton instance pointer for timing access.
        static AVBDPhysicsServer3D *singleton;
        bool debug_contacts_requested = false;
        int32_t debug_contacts_max = 0;

        // Cached direct space state handed to the engine; owned, freed with the space.
        AVBDDirectSpaceState3D *direct_state = nullptr;

        // Solver-body -> server id, rebuilt with the space. The direct space state's
        // queries get their candidates back as Rigid* from the BVH; this is the O(1)
        // bridge back to RIDs (body_id_of used to scan every body per result).
        ankerl::unordered_dense::map<const avbd::Rigid *, uint64_t> rigid_ids;
    };

    // One area, as the engine described it. Monitored each step against the bodies of
    // its space (see _monitor_areas).
    struct AreaData {
        RID space;
        std::vector<ShapeData> shapes;
        std::vector<Transform3D> shape_transforms;
        std::vector<bool> shape_disabled;
        Transform3D transform;
        uint64_t instance_id = 0;
        uint32_t collision_layer = 1;
        uint32_t collision_mask = 1;
        bool monitorable = false;
        bool monitor_callback_set = false;
        bool area_monitor_callback_set = false;
        bool ray_pickable = true;
        Callable monitor_callback;        // bodies entering/leaving this area
        Callable area_monitor_callback;   // other areas entering/leaving this area
        // Parameters recorded for the getter; they do not act on the solver this round.
        ankerl::unordered_dense::map<int, Variant> params;
        // Solver ids inside the area after the last step, for enter/exit diffing.
        std::vector<uint64_t> prev_inside_bodies;
        std::vector<uint64_t> prev_inside_areas;
    };

    // Object kinds are kept in separate maps so an id that arrives for the wrong kind is a miss
    // rather than a reinterpretation. `ankerl::unordered_dense` stores values contiguously and
    // relocates them on rehash, which is fine for the plain-data records. SpaceData owns a
    // Solver - a non-movable RAII holder of intrusive lists whose nodes point back at it - so
    // spaces keeps it behind a stable unique_ptr: the flat map only ever moves the handle,
    // never the solver, and the space's address survives any rehash.
    ankerl::unordered_dense::map<uint64_t, std::unique_ptr<SpaceData>> spaces;
    ankerl::unordered_dense::map<uint64_t, BodyData> bodies;
    ankerl::unordered_dense::map<uint64_t, ShapeData> shapes;
    ankerl::unordered_dense::map<uint64_t, JointData> joints;
    ankerl::unordered_dense::map<uint64_t, AreaData> areas;

    uint64_t next_id = 1;

    // Godot pauses the SceneTree with PhysicsServer3D.set_active(false); stepping must
    // follow the game, so _step refuses to run while the engine says the server is inactive.
    bool server_active = true;
    bool flushing_queries = false;
    // Set once any space has stepped. _flush_queries delivers state to the scene only
    // after that, so the engine never reads a body pose the solver has not produced.
    bool has_stepped = false;

    // Deepest contact of the last completed step, snapshotted on the stepping thread.
    // The getters run on the main thread, where walking solver.forces races with the
    // solver's manifold creation/revival/reaping and can report depths for pooled
    // manifolds that no longer describe a touching pair.
    double diag_max_penetration = 0.0;
    uint64_t diag_deepest_a = 0;
    uint64_t diag_deepest_b = 0;

    // Cached state objects handed to the engine's callbacks; re-aimed per body before use.
    AVBDDirectBodyState3D *body_state = nullptr;

    uint64_t id_of(const RID &p_rid) const;
    RID make_rid(uint64_t p_id) const;

    BodyData *find_body(const RID &p_rid);
    const BodyData *find_body(const RID &p_rid) const;
    ShapeData *find_shape(const RID &p_rid);
    SpaceData *find_space(const RID &p_rid);
    AreaData *find_area(const RID &p_rid);

    // Godot keeps gravity in the project settings rather than handing it to the server.
    static void read_project_gravity(SpaceData &p_space);
    // Solver parameters live in project settings the same way ("physics/avbd/*").
    void read_project_params(SpaceData &p_space);
    // ProjectSettings::settings_changed lands here (deferred); marks every space's
    // parameters for a re-read on the next step.
    void _on_project_settings_changed();

    // Rebuild the solver for one space from the bodies and joints assigned to it. Needed
    // whenever a body, its shapes, its mode or a joint changes, because those decide what the
    // solver holds.
    // The source location is the rebuild reason: every call site is a distinct engine
    // entry point, and knowing which one fired (and when) is what tells a mid-run
    // rebuild - which destroys every manifold - apart from scene setup.
    void rebuild_space(const RID &p_space,
            std::source_location where = std::source_location::current());

    // Spaces whose body set changed since the last step; rebuilt at the top of _step so
    // a body joining mid-setup lands with its shapes attached.
    ankerl::unordered_dense::map<uint64_t, bool> space_rebuild_pending;
    // Lifetime count of rebuild_space() calls; reported by get_solver_sleep_stats().
    uint64_t rebuild_total = 0;
    // Why the rebuilds happened, and when the last one landed relative to _step() calls. A
    // rebuild destroys every manifold, so one arriving mid-run silently discards all
    // warm-started contact state; these make that visible from a script.
    ankerl::unordered_dense::map<std::string, uint64_t> rebuild_reasons;
    uint64_t step_total = 0;
    uint64_t rebuild_last_step = 0;
    String rebuild_last_reason;

    // Density the solver should see for a body: mass / volume of the first shape, or the house
    // density when no mass was set.
    static float density_for(const BodyData &p_body);

    // Joint mappings: build one solver GenericJoint per engine joint record.
    void rebuild_6dof_joint(SpaceData *p_space, const JointData &p_joint, avbd::Rigid *p_rigid_a, avbd::Rigid *p_rigid_b);
    void rebuild_pin_joint(SpaceData *p_space, const JointData &p_joint, avbd::Rigid *p_rigid_a, avbd::Rigid *p_rigid_b);
    void rebuild_hinge_joint(SpaceData *p_space, const JointData &p_joint, avbd::Rigid *p_rigid_a, avbd::Rigid *p_rigid_b);
    void rebuild_slider_joint(SpaceData *p_space, const JointData &p_joint, avbd::Rigid *p_rigid_a, avbd::Rigid *p_rigid_b);
    void rebuild_cone_twist_joint(SpaceData *p_space, const JointData &p_joint, avbd::Rigid *p_rigid_a, avbd::Rigid *p_rigid_b);
    void _rebuild_for_joint(const RID &p_joint);

    // --- DirectBodyState accessors (used by AVBDDirectBodyState3D) ---------------
    // On the server rather than the state object because they need the body map.
public:
    // --- query helpers for AVBDDirectSpaceState3D -------------------------------
    // One simulated body a query may hit.
    struct QueryCandidate {
        avbd::Rigid *rigid;
        uint64_t body_id; // server map id
    };
    // Every simulated body of the space whose layer survives the query mask.
    std::vector<QueryCandidate> query_candidates(const RID &p_space, uint32_t p_mask) const;
    // The same contract, narrowed by the last built broad-phase tree: only bodies whose
    // bound overlaps the query AABB (sim space). Before any step has built a tree this
    // answers from the body list with per-body bounds.
    std::vector<QueryCandidate> query_candidates_aabb(const RID &p_space, avbd::float3 p_min,
            avbd::float3 p_max, uint32_t p_mask) const;

    // Solver::pick with a mask, on one space's solver.
    avbd::Rigid *solver_pick(const RID &p_space, avbd::float3 p_origin, avbd::float3 p_dir, uint32_t p_mask,
            avbd::float3 &r_local);
    // The body record id owning a solver rigid within a space (0 if none).
    uint64_t body_id_of(const RID &p_space, const avbd::Rigid *p_rigid) const;
    RID solver_pick_rid(const RID &p_space, const avbd::Rigid *p_rigid);
    uint64_t body_instance_id(const RID &p_space, const avbd::Rigid *p_rigid);
    Vector3 body_velocity_of(const RID &p_space, const avbd::Rigid *p_rigid);
    const ShapeData *shape_data(const RID &p_shape) const;

    void body_state_add_constant_force(const RID &p_body, const Vector3 &p_force, const Vector3 &p_position);
    void body_state_add_constant_torque(const RID &p_body, const Vector3 &p_torque);
    void body_state_set_constant_force(const RID &p_body, const Vector3 &p_force);
    Vector3 body_state_get_constant_force(const RID &p_body) const;
    void body_state_set_constant_torque(const RID &p_body, const Vector3 &p_torque);
    Vector3 body_state_get_constant_torque(const RID &p_body) const;
    int32_t body_state_contact_count(const RID &p_body) const;
    Vector3 body_state_contact_position(const RID &p_body, int32_t p_index) const;
    RID body_state_contact_collider(const RID &p_body, int32_t p_index) const;
    uint64_t body_state_contact_collider_id(const RID &p_body, int32_t p_index) const;
    Vector3 body_state_contact_collider_velocity(const RID &p_body, int32_t p_index) const;
    double body_state_step() const;
    PhysicsDirectSpaceState3D *body_state_space_state(const RID &p_body);

    // Solver performance timing (ms).
    double get_solver_broad_phase_time() const;
    double get_solver_colour_graph_time() const;
    double get_solver_solve_time() const;
    double get_solver_finish_time() const;

    // Penetration diagnostics for the headless stress suites. A manifold's contact error
    // C0 carries the geometric overlap captured at the start of the step: its normal
    // component is `COLLISION_MARGIN - depth`, so the depth in metres is margin - C0.x.
    // Returns the deepest contact across every space (0 when nothing penetrates).
    double get_solver_max_penetration() const;
    // The deepest contact as [body_a_instance_id, body_b_instance_id, depth] (Godot
    // instance ids, so a script can resolve the nodes). Empty when nothing penetrates.
    Array get_solver_deepest_contact() const;
    // Sleep/wake diagnostics for the stress suites. A sleeping body is skipped by the
    // velocity update in Solver::finishBody, and a manifold with neither side awake is
    // frozen (its contact state is deliberately left stale), so both counts say whether a
    // penetration measurement was taken against a live solve or a parked one. `rebuilds` is
    // cumulative over the server's lifetime: a rebuild destroys every manifold, so a rising
    // count during a run means warm-started contact state is being thrown away.
    Dictionary get_solver_sleep_stats() const;
    // Get the singleton instance of AVBDPhysicsServer3D.
    static AVBDPhysicsServer3D *get_singleton();

protected:
    static void _bind_methods() {
        // Solver performance timing
        ClassDB::bind_method(D_METHOD("get_solver_broad_phase_time"), &AVBDPhysicsServer3D::get_solver_broad_phase_time);
        ClassDB::bind_method(D_METHOD("get_solver_colour_graph_time"), &AVBDPhysicsServer3D::get_solver_colour_graph_time);
        ClassDB::bind_method(D_METHOD("get_solver_solve_time"), &AVBDPhysicsServer3D::get_solver_solve_time);
        ClassDB::bind_method(D_METHOD("get_solver_finish_time"), &AVBDPhysicsServer3D::get_solver_finish_time);
        ClassDB::bind_method(D_METHOD("get_solver_max_penetration"), &AVBDPhysicsServer3D::get_solver_max_penetration);
        ClassDB::bind_method(D_METHOD("get_solver_deepest_contact"), &AVBDPhysicsServer3D::get_solver_deepest_contact);
        ClassDB::bind_method(D_METHOD("get_solver_sleep_stats"), &AVBDPhysicsServer3D::get_solver_sleep_stats);
        // Target of the ProjectSettings::settings_changed connection made in _init. Bound
        // so the StringName callable can reach it; not meant for scripts.
        ClassDB::bind_method(D_METHOD("_on_project_settings_changed"), &AVBDPhysicsServer3D::_on_project_settings_changed);
        // ClassDB::bind_method(D_METHOD("get_singleton"), &AVBDPhysicsServer3D::get_singleton); // static method, not bound
    }

public:
    AVBDPhysicsServer3D() = default;
    ~AVBDPhysicsServer3D() override;

    // --- lifecycle ---------------------------------------------------------------
    void _init() override;
    void _finish() override;
    void _step(float p_step) override;
    void _sync() override;
    void _end_sync() override;
    void _flush_queries() override;
    void _free_rid(const RID &p_rid) override;
    void _set_active(bool p_active) override;
    bool _is_flushing_queries() const override;
    int32_t _get_process_info(PhysicsServer3D::ProcessInfo p_info) override;

    // --- spaces ------------------------------------------------------------------
    RID _space_create() override;
    void _space_set_active(const RID &p_space, bool p_active) override;
    bool _space_is_active(const RID &p_space) const override;
    void _space_set_param(const RID &p_space, PhysicsServer3D::SpaceParameter p_param, float p_value) override;
    float _space_get_param(const RID &p_space, PhysicsServer3D::SpaceParameter p_param) const override;
    PhysicsDirectSpaceState3D *_space_get_direct_state(const RID &p_space) override;
    void _space_set_debug_contacts(const RID &p_space, int32_t p_max_contacts) override;
    PackedVector3Array _space_get_contacts(const RID &p_space) const override;
    int32_t _space_get_contact_count(const RID &p_space) const override;

    // --- shapes ------------------------------------------------------------------
        void _shape_set_custom_solver_bias(const RID &p_shape, float p_bias) override;
    float _shape_get_custom_solver_bias(const RID &p_shape) const override;
    void _shape_set_margin(const RID &p_shape, float p_margin) override;
    float _shape_get_margin(const RID &p_shape) const override;
    RID _capsule_shape_create() override;
    RID _convex_polygon_shape_create() override;
    RID _concave_polygon_shape_create() override;
    RID _heightmap_shape_create() override;
    RID _world_boundary_shape_create() override;
    RID _separation_ray_shape_create() override;
    RID _custom_shape_create() override;
    bool _body_test_motion(const RID &p_body, const Transform3D &p_from, const Vector3 &p_motion, float p_margin,
            int32_t p_max_collisions, bool p_collide_separation_ray, bool p_recovery_as_collision,
            PhysicsServer3DExtensionMotionResult *p_result) const override;
RID _box_shape_create() override;
    RID _sphere_shape_create() override;
    RID _cylinder_shape_create() override;
    void _shape_set_data(const RID &p_shape, const Variant &p_data) override;
    PhysicsServer3D::ShapeType _shape_get_type(const RID &p_shape) const override;
    Variant _shape_get_data(const RID &p_shape) const override;

    // --- bodies ------------------------------------------------------------------
    RID _body_create() override;
    void _body_set_space(const RID &p_body, const RID &p_space) override;
    RID _body_get_space(const RID &p_body) const override;
    void _body_set_mode(const RID &p_body, PhysicsServer3D::BodyMode p_mode) override;
    PhysicsServer3D::BodyMode _body_get_mode(const RID &p_body) const override;
    void _body_set_state(const RID &p_body, PhysicsServer3D::BodyState p_state, const Variant &p_value) override;
    Variant _body_get_state(const RID &p_body, PhysicsServer3D::BodyState p_state) const override;
    void _body_set_state_sync_callback(const RID &p_body, const Callable &p_callable) override;
    void _body_attach_object_instance_id(const RID &p_body, uint64_t p_id) override;
    uint64_t _body_get_object_instance_id(const RID &p_body) const override;

    // --- soft-body server API: documented unsupported this round -----------------
    // The engine calls these only for SoftBody3D nodes, which the demo scenes do not
    // use (the AVBDSoftBody3D node class runs its own solver). Neutral answers keep a
    // stray SoftBody3D from erroring; nothing is simulated.
    RID _soft_body_create() override;
    void _soft_body_set_space(const RID &p_body, const RID &p_space) override;
    RID _soft_body_get_space(const RID &p_body) const override;
    void _soft_body_set_mesh(const RID &p_body, const RID &p_mesh) override;
    void _soft_body_set_transform(const RID &p_body, const Transform3D &p_transform) override;
    void _soft_body_set_collision_layer(const RID &p_body, uint32_t p_layer) override;
    uint32_t _soft_body_get_collision_layer(const RID &p_body) const override;
    void _soft_body_set_collision_mask(const RID &p_body, uint32_t p_mask) override;
    uint32_t _soft_body_get_collision_mask(const RID &p_body) const override;
    void _soft_body_add_collision_exception(const RID &p_body, const RID &p_body_b) override;
    void _soft_body_remove_collision_exception(const RID &p_body, const RID &p_body_b) override;
    TypedArray<RID> _soft_body_get_collision_exceptions(const RID &p_body) const override;
    void _soft_body_set_total_mass(const RID &p_body, float p_total_mass) override;
    float _soft_body_get_total_mass(const RID &p_body) const override;
    void _soft_body_set_linear_stiffness(const RID &p_body, float p_linear_stiffness) override;
    float _soft_body_get_linear_stiffness(const RID &p_body) const override;
    void _soft_body_set_pressure_coefficient(const RID &p_body, float p_pressure_coefficient) override;
    float _soft_body_get_pressure_coefficient(const RID &p_body) const override;
    void _soft_body_set_damping_coefficient(const RID &p_body, float p_damping_coefficient) override;
    float _soft_body_get_damping_coefficient(const RID &p_body) const override;
    void _soft_body_set_drag_coefficient(const RID &p_body, float p_drag_coefficient) override;
    float _soft_body_get_drag_coefficient(const RID &p_body) const override;
    void _soft_body_set_simulation_precision(const RID &p_body, int32_t p_simulation_precision) override;
    int32_t _soft_body_get_simulation_precision(const RID &p_body) const override;
    void _soft_body_set_shrinking_factor(const RID &p_body, float p_shrinking_factor) override;
    float _soft_body_get_shrinking_factor(const RID &p_body) const override;
    void _soft_body_set_state(const RID &p_body, PhysicsServer3D::BodyState p_state, const Variant &p_variant) override;
    Variant _soft_body_get_state(const RID &p_body, PhysicsServer3D::BodyState p_state) const override;
    void _soft_body_apply_central_force(const RID &p_body, const Vector3 &p_force) override;
    void _soft_body_apply_central_impulse(const RID &p_body, const Vector3 &p_impulse) override;
    void _soft_body_apply_point_force(const RID &p_body, int32_t p_point_index, const Vector3 &p_force) override;
    void _soft_body_apply_point_impulse(const RID &p_body, int32_t p_point_index, const Vector3 &p_impulse) override;
    void _soft_body_move_point(const RID &p_body, int32_t p_point_index, const Vector3 &p_global_position) override;
    void _soft_body_pin_point(const RID &p_body, int32_t p_point_index, bool p_pin) override;
    bool _soft_body_is_point_pinned(const RID &p_body, int32_t p_point_index) const override;
    Vector3 _soft_body_get_point_global_position(const RID &p_body, int32_t p_point_index) const override;
    void _soft_body_remove_all_pinned_points(const RID &p_body) override;
    void _soft_body_set_ray_pickable(const RID &p_body, bool p_enable) override;
    AABB _soft_body_get_bounds(const RID &p_body) const override;
    void _soft_body_update_rendering_server(const RID &p_body,
            PhysicsServer3DRenderingServerHandler *p_rendering_server_handler) override;
    void _body_set_ray_pickable(const RID &p_body, bool p_enable) override;
    bool _body_is_continuous_collision_detection_enabled(const RID &p_body) const override;
    void _body_set_enable_continuous_collision_detection(const RID &p_body, bool p_enable) override;
    void _body_set_collision_layer(const RID &p_body, uint32_t p_layer) override;
    uint32_t _body_get_collision_layer(const RID &p_body) const override;
    void _body_set_collision_mask(const RID &p_body, uint32_t p_mask) override;
    uint32_t _body_get_collision_mask(const RID &p_body) const override;
    void _body_set_collision_priority(const RID &p_body, float p_priority) override;
    float _body_get_collision_priority(const RID &p_body) const override;
    void _body_set_user_flags(const RID &p_body, uint32_t p_flags) override;
    uint32_t _body_get_user_flags(const RID &p_body) const override;
    void _body_reset_mass_properties(const RID &p_body) override;
    void _body_add_collision_exception(const RID &p_body, const RID &p_excepted_body) override;
    void _body_remove_collision_exception(const RID &p_body, const RID &p_excepted_body) override;
    TypedArray<RID> _body_get_collision_exceptions(const RID &p_body) const override;
    void _body_set_max_contacts_reported(const RID &p_body, int32_t p_amount) override;
    int32_t _body_get_max_contacts_reported(const RID &p_body) const override;
    void _body_set_contacts_reported_depth_threshold(const RID &p_body, float p_threshold) override;
    float _body_get_contacts_reported_depth_threshold(const RID &p_body) const override;
    void _body_set_omit_force_integration(const RID &p_body, bool p_enable) override;
    bool _body_is_omitting_force_integration(const RID &p_body) const override;
    void _body_set_force_integration_callback(const RID &p_body, const Callable &p_callable, const Variant &p_userdata) override;
    void _body_set_constant_force(const RID &p_body, const Vector3 &p_force) override;
    Vector3 _body_get_constant_force(const RID &p_body) const override;
    void _body_set_constant_torque(const RID &p_body, const Vector3 &p_torque) override;
    Vector3 _body_get_constant_torque(const RID &p_body) const override;
    void _body_set_axis_velocity(const RID &p_body, const Vector3 &p_axis_velocity) override;
    bool _body_is_axis_locked(const RID &p_body, PhysicsServer3D::BodyAxis p_axis) const override;
    PhysicsDirectBodyState3D *_body_get_direct_state(const RID &p_body) override;

    // --- body parameters and forces ----------------------------------------------
    void _body_set_param(const RID &p_body, PhysicsServer3D::BodyParameter p_param, const Variant &p_value) override;
    Variant _body_get_param(const RID &p_body, PhysicsServer3D::BodyParameter p_param) const override;
    void _body_apply_central_impulse(const RID &p_body, const Vector3 &p_impulse) override;
    void _body_apply_impulse(const RID &p_body, const Vector3 &p_impulse, const Vector3 &p_position) override;
    void _body_apply_torque_impulse(const RID &p_body, const Vector3 &p_impulse) override;
    void _body_apply_central_force(const RID &p_body, const Vector3 &p_force) override;
    void _body_apply_force(const RID &p_body, const Vector3 &p_force, const Vector3 &p_position) override;
    void _body_apply_torque(const RID &p_body, const Vector3 &p_torque) override;
    void _body_add_constant_central_force(const RID &p_body, const Vector3 &p_force) override;
    void _body_add_constant_force(const RID &p_body, const Vector3 &p_force, const Vector3 &p_position) override;
    void _body_add_constant_torque(const RID &p_body, const Vector3 &p_torque) override;
    void _body_set_axis_lock(const RID &p_body, PhysicsServer3D::BodyAxis p_axis, bool p_lock) override;

    void _body_add_shape(const RID &p_body, const RID &p_shape, const Transform3D &p_transform,
            bool p_disabled) override;
    void _body_set_shape(const RID &p_body, int32_t p_shape_idx, const RID &p_shape) override;
    void _body_set_shape_transform(const RID &p_body, int32_t p_shape_idx, const Transform3D &p_transform) override;
    void _body_set_shape_disabled(const RID &p_body, int32_t p_shape_idx, bool p_disabled) override;
    int32_t _body_get_shape_count(const RID &p_body) const override;
    RID _body_get_shape(const RID &p_body, int32_t p_shape_idx) const override;
    Transform3D _body_get_shape_transform(const RID &p_body, int32_t p_shape_idx) const override;
    void _body_remove_shape(const RID &p_body, int32_t p_shape_idx) override;
    void _body_clear_shapes(const RID &p_body) override;

    // --- joints -------------------------------------------------------------------
    // Only the generic 6-DOF joint gets a solver counterpart: AVBD's GenericJoint covers every
    // configuration the car and the demos need, and the other joint types are accepted and
    // tracked (so scenes load) without acting. Their setters record and their getters replay.
    RID _joint_create() override;
    void _joint_clear(const RID &p_joint) override;
    PhysicsServer3D::JointType _joint_get_type(const RID &p_joint) const override;
    void _joint_set_solver_priority(const RID &p_joint, int32_t p_priority) override;
    int32_t _joint_get_solver_priority(const RID &p_joint) const override;
    void _joint_disable_collisions_between_bodies(const RID &p_joint, bool p_disable) override;
    bool _joint_is_disabled_collisions_between_bodies(const RID &p_joint) const override;

    void _joint_make_generic_6dof(const RID &p_joint, const RID &p_body_a, const Transform3D &p_local_ref_a,
            const RID &p_body_b, const Transform3D &p_local_ref_b) override;
    void _generic_6dof_joint_set_param(const RID &p_joint, Vector3::Axis p_axis,
            PhysicsServer3D::G6DOFJointAxisParam p_param, float p_value) override;
    float _generic_6dof_joint_get_param(const RID &p_joint, Vector3::Axis p_axis,
            PhysicsServer3D::G6DOFJointAxisParam p_param) const override;
    void _generic_6dof_joint_set_flag(const RID &p_joint, Vector3::Axis p_axis,
            PhysicsServer3D::G6DOFJointAxisFlag p_flag, bool p_enable) override;
    bool _generic_6dof_joint_get_flag(const RID &p_joint, Vector3::Axis p_axis,
            PhysicsServer3D::G6DOFJointAxisFlag p_flag) const override;

    // --- joints: the other four types, mapped onto avbd::GenericJoint -------------
    void _joint_make_pin(const RID &p_joint, const RID &p_body_a, const Vector3 &p_local_a,
            const RID &p_body_b, const Vector3 &p_local_b) override;
    void _pin_joint_set_param(const RID &p_joint, PhysicsServer3D::PinJointParam p_param, float p_value) override;
    float _pin_joint_get_param(const RID &p_joint, PhysicsServer3D::PinJointParam p_param) const override;
    void _pin_joint_set_local_a(const RID &p_joint, const Vector3 &p_local_a) override;
    Vector3 _pin_joint_get_local_a(const RID &p_joint) const override;
    void _pin_joint_set_local_b(const RID &p_joint, const Vector3 &p_local_b) override;
    Vector3 _pin_joint_get_local_b(const RID &p_joint) const override;

    void _joint_make_hinge(const RID &p_joint, const RID &p_body_a, const Transform3D &p_hinge_a,
            const RID &p_body_b, const Transform3D &p_hinge_b) override;
    void _joint_make_hinge_simple(const RID &p_joint, const RID &p_body_a, const Vector3 &p_pivot_a,
            const Vector3 &p_axis_a, const RID &p_body_b, const Vector3 &p_pivot_b,
            const Vector3 &p_axis_b) override;
    void _hinge_joint_set_param(const RID &p_joint, PhysicsServer3D::HingeJointParam p_param, float p_value) override;
    float _hinge_joint_get_param(const RID &p_joint, PhysicsServer3D::HingeJointParam p_param) const override;
    void _hinge_joint_set_flag(const RID &p_joint, PhysicsServer3D::HingeJointFlag p_flag, bool p_enabled) override;
    bool _hinge_joint_get_flag(const RID &p_joint, PhysicsServer3D::HingeJointFlag p_flag) const override;

    void _joint_make_slider(const RID &p_joint, const RID &p_body_a, const Transform3D &p_local_ref_a,
            const RID &p_body_b, const Transform3D &p_local_ref_b) override;
    void _slider_joint_set_param(const RID &p_joint, PhysicsServer3D::SliderJointParam p_param, float p_value) override;
    float _slider_joint_get_param(const RID &p_joint, PhysicsServer3D::SliderJointParam p_param) const override;

    void _joint_make_cone_twist(const RID &p_joint, const RID &p_body_a, const Transform3D &p_local_ref_a,
            const RID &p_body_b, const Transform3D &p_local_ref_b) override;
    void _cone_twist_joint_set_param(const RID &p_joint, PhysicsServer3D::ConeTwistJointParam p_param,
            float p_value) override;
    float _cone_twist_joint_get_param(const RID &p_joint, PhysicsServer3D::ConeTwistJointParam p_param) const override;

    // --- areas --------------------------------------------------------------------
    // Created and tracked; every storage virtual is implemented and monitoring runs each
    // step (see _monitor_areas).
    RID _area_create() override;
    void _area_set_space(const RID &p_area, const RID &p_space) override;
    RID _area_get_space(const RID &p_area) const override;
    void _area_set_transform(const RID &p_area, const Transform3D &p_transform) override;
    Transform3D _area_get_transform(const RID &p_area) const override;
    void _area_add_shape(const RID &p_area, const RID &p_shape, const Transform3D &p_transform, bool p_disabled)
            override;
    void _area_set_shape(const RID &p_area, int32_t p_shape_idx, const RID &p_shape) override;
    void _area_set_shape_transform(const RID &p_area, int32_t p_shape_idx, const Transform3D &p_transform) override;
    void _area_set_shape_disabled(const RID &p_area, int32_t p_shape_idx, bool p_disabled) override;
    int32_t _area_get_shape_count(const RID &p_area) const override;
    RID _area_get_shape(const RID &p_area, int32_t p_shape_idx) const override;
    Transform3D _area_get_shape_transform(const RID &p_area, int32_t p_shape_idx) const override;
    void _area_remove_shape(const RID &p_area, int32_t p_shape_idx) override;
    void _area_clear_shapes(const RID &p_area) override;
    void _area_set_param(const RID &p_area, PhysicsServer3D::AreaParameter p_param, const Variant &p_value) override;
    Variant _area_get_param(const RID &p_area, PhysicsServer3D::AreaParameter p_param) const override;
    void _area_set_collision_layer(const RID &p_area, uint32_t p_layer) override;
    uint32_t _area_get_collision_layer(const RID &p_area) const override;
    void _area_set_collision_mask(const RID &p_area, uint32_t p_mask) override;
    uint32_t _area_get_collision_mask(const RID &p_area) const override;
    void _area_set_monitorable(const RID &p_area, bool p_monitorable) override;
    void _area_set_ray_pickable(const RID &p_area, bool p_enable) override;
    void _area_attach_object_instance_id(const RID &p_area, uint64_t p_id) override;
    uint64_t _area_get_object_instance_id(const RID &p_area) const override;
    void _area_set_monitor_callback(const RID &p_area, const Callable &p_callback) override;
    void _area_set_area_monitor_callback(const RID &p_area, const Callable &p_callback) override;

    // --- per-step area monitoring ---------------------------------------------------
    // Diff every area against the bodies/areas of its space and fire the callbacks.
    void _monitor_areas();
};

} // namespace godot

#endif // AVBD_PHYSICS_SERVER3D_HPP

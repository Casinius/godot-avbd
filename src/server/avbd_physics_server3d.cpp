/*
 * AVBD as a Godot physics server. See the header for the design and the unsupported list.
 */

#include "server/avbd_physics_server3d.hpp"

#include "nodes/godot_convert.hpp"
#include "server/avbd_direct_body_state3d.hpp"
#include "server/avbd_direct_space_state3d.hpp"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <tuple>

#include <godot_cpp/classes/box_shape3d.hpp>
#include <godot_cpp/classes/cylinder_shape3d.hpp>
#include <godot_cpp/classes/sphere_shape3d.hpp>
#include <godot_cpp/classes/project_settings.hpp>
#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/variant/array.hpp>

using namespace godot;
using namespace avbd_godot;

// Godot 4.3 added per-axis springs to the generic 6-DOF joint and grew the server enums; the
// godot-cpp headers this extension builds against predate that, so the new slots are named
// here with the engine's values (servers/physics_server_3d.h, G6DOFJointAxisParam/Flag).
namespace g6dof_spring {
constexpr int LINEAR_STIFFNESS = 7;
constexpr int LINEAR_DAMPING = 8;
constexpr int LINEAR_EQUILIBRIUM = 9;
constexpr int ANGULAR_STIFFNESS = 19;
constexpr int ANGULAR_DAMPING = 20;
constexpr int ANGULAR_EQUILIBRIUM = 21;
constexpr int FLAG_ENABLE_ANGULAR = 2;
constexpr int FLAG_ENABLE_LINEAR = 3;
} // namespace g6dof_spring

// -----------------------------------------------------------------------------
// Handles
//
// The extension interface offers no way to allocate a RID, so the server mints its own: an
// integer id written into the RID's opaque bytes through the pointer godot-cpp does expose. The
// engine only ever hands these values back to us, so the encoding is ours to choose - but it must
// round-trip, which is what id_of() checks by reading them back.
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::make_rid(uint64_t p_id) const {
    RID rid;
    std::memcpy(rid._native_ptr(), &p_id, sizeof(p_id));
    return rid;
}

uint64_t AVBDPhysicsServer3D::id_of(const RID &p_rid) const {
    uint64_t id = 0;
    std::memcpy(&id, p_rid._native_ptr(), sizeof(id));
    return id;
}

AVBDPhysicsServer3D::BodyData *AVBDPhysicsServer3D::find_body(const RID &p_rid) {
    const auto found = bodies.find(id_of(p_rid));
    return found == bodies.end() ? nullptr : &found->second;
}

const AVBDPhysicsServer3D::BodyData *AVBDPhysicsServer3D::find_body(const RID &p_rid) const {
    const auto found = bodies.find(id_of(p_rid));
    return found == bodies.end() ? nullptr : &found->second;
}

AVBDPhysicsServer3D::ShapeData *AVBDPhysicsServer3D::find_shape(const RID &p_rid) {
    const auto found = shapes.find(id_of(p_rid));
    return found == shapes.end() ? nullptr : &found->second;
}

AVBDPhysicsServer3D::SpaceData *AVBDPhysicsServer3D::find_space(const RID &p_rid) {
    const auto found = spaces.find(id_of(p_rid));
    return found == spaces.end() ? nullptr : &found->second;
}

// -----------------------------------------------------------------------------
// Lifecycle
// -----------------------------------------------------------------------------

// Spaces may outlive _finish in edge paths (a space freed without the engine stopping
// first): the destructor is the last chance to release the cached state objects.
AVBDPhysicsServer3D::~AVBDPhysicsServer3D() {
    for (auto &[id, space] : spaces) {
        (void)id;
        if (space.direct_state != nullptr) {
            memdelete(space.direct_state);
            space.direct_state = nullptr;
        }
    }
    if (body_state != nullptr) {
        memdelete(body_state);
        body_state = nullptr;
    }
}

void AVBDPhysicsServer3D::_finish() {
    // The engine is shutting the server down; drop everything it owns. The cached state
    // objects are raw pointers inside the records, so they go first.
    for (auto &[id, space] : spaces) {
        (void)id;
        if (space.direct_state != nullptr) {
            memdelete(space.direct_state);
            space.direct_state = nullptr;
        }
    }
    spaces.clear();
    bodies.clear();
    shapes.clear();
    joints.clear();
    areas.clear();
    if (body_state != nullptr) {
        memdelete(body_state);
        body_state = nullptr;
    }
}

void AVBDPhysicsServer3D::_init() {
    // Spaces created after this moment read gravity themselves; this covers the ones the engine
    // creates before any of our code runs.
    for (auto &[id, space] : spaces) {
        read_project_gravity(space);
        read_project_params(space);
    }
}

void AVBDPhysicsServer3D::_flush_queries() {
    flushing_queries = true;
    // Nothing accumulates queries yet; the direct space state answers immediately.
    flushing_queries = false;
}

bool AVBDPhysicsServer3D::_is_flushing_queries() const {
    return flushing_queries;
}

int32_t AVBDPhysicsServer3D::_get_process_info(PhysicsServer3D::ProcessInfo p_info) {
    switch (p_info) {
        case PhysicsServer3D::INFO_ACTIVE_OBJECTS: {
            int32_t count = 0;
            for (const auto &[id, body] : bodies) {
                (void)id;
                if (body.rigid != nullptr) {
                    count++;
                }
            }
            return count;
        }
        case PhysicsServer3D::INFO_ISLAND_COUNT:
            return static_cast<int32_t>(spaces.size());
        case PhysicsServer3D::INFO_COLLISION_PAIRS: {
            // Manifolds are the forces with contact points.
            int32_t count = 0;
            for (const auto &[id, space] : spaces) {
                (void)id;
                for (const avbd::Force *force = space.solver.forces; force != nullptr; force = force->next) {
                    if (force->contactPointCount() > 0) {
                        count++;
                    }
                }
            }
            return count;
        }
        default:
            return 0;
    }
}

void AVBDPhysicsServer3D::_set_active(bool p_active) {
    // SceneTree::set_pause calls PhysicsServer3D.set_active(!paused): this is the point
    // where simulation must follow the game. _step refuses to run while false.
    server_active = p_active;
}

void AVBDPhysicsServer3D::_sync() {
    // Called from the physics thread before stepping. The bodies were already written by
    // _body_set_state, so there is nothing to pull in here.
}

void AVBDPhysicsServer3D::_end_sync() {
    // ...and after stepping, which is when a real implementation would push results back. Ours
    // pushes them at the end of _step, where the solver's output is in hand.
}

void AVBDPhysicsServer3D::_free_rid(const RID &p_rid) {
    const uint64_t id = id_of(p_rid);

    const auto space_found = spaces.find(id);
    if (space_found != spaces.end()) {
        if (space_found->second.direct_state != nullptr) {
            memdelete(space_found->second.direct_state);
        }
        spaces.erase(space_found);
        return;
    }
    if (bodies.erase(id) > 0) {
        return;
    }
    if (shapes.erase(id) > 0) {
        return;
    }
    if (joints.erase(id) > 0) {
        return;
    }
    areas.erase(id);
}

// -----------------------------------------------------------------------------
// Spaces
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_space_create() {
    const uint64_t id = next_id++;
    // emplace, not assignment: Solver owns intrusive lists and is deliberately non-movable.
    SpaceData &space = spaces.emplace(std::piecewise_construct, std::forward_as_tuple(id),
            std::forward_as_tuple())
                               .first->second;
    read_project_gravity(space);
    read_project_params(space);
    return make_rid(id);
}

// Solver parameters come from project settings the same way gravity does, so a project can
// tune the solver without a node layer. Keys absent from a project fall back to the solver's
// defaults.
void AVBDPhysicsServer3D::read_project_params(SpaceData &p_space) {
    const ProjectSettings *settings = ProjectSettings::get_singleton();
    if (settings == nullptr) {
        return;
    }

    p_space.threads = static_cast<int>(settings->get_setting("physics/avbd/threads", 0));
    p_space.iterations = static_cast<int>(settings->get_setting("physics/avbd/iterations", 10));
    p_space.alpha = settings->get_setting("physics/avbd/alpha", 0.99);
    p_space.beta_linear = settings->get_setting("physics/avbd/beta_linear", 10000.0);
    p_space.beta_angular = settings->get_setting("physics/avbd/beta_angular", 100.0);
    p_space.gamma = settings->get_setting("physics/avbd/gamma", 0.999);
    p_space.substeps = static_cast<int>(settings->get_setting("physics/avbd/substeps", 1));
}

void AVBDPhysicsServer3D::_space_set_active(const RID &p_space, bool p_active) {
    if (SpaceData *space = find_space(p_space)) {
        space->active = p_active;
    }
}

bool AVBDPhysicsServer3D::_space_is_active(const RID &p_space) const {
    const auto found = spaces.find(id_of(p_space));
    return found == spaces.end() ? false : found->second.active;
}

void AVBDPhysicsServer3D::_space_set_param(const RID &p_space, PhysicsServer3D::SpaceParameter p_param, double p_value) {
    // The space parameters are contact tolerances and sleep thresholds, none of which map onto
    // AVBD's model. They are accepted and ignored rather than rejected, so a project's settings do
    // not error out on load. Gravity does not come through here at all: Godot keeps it in the
    // project settings, and a server is expected to read it (see read_project_gravity).
    (void)p_space;
    (void)p_param;
    (void)p_value;
}

// Godot stores gravity in the project settings rather than passing it to the server, so the server
// reads it - the same thing the built-in one does.
void AVBDPhysicsServer3D::read_project_gravity(SpaceData &p_space) {
    const ProjectSettings *settings = ProjectSettings::get_singleton();
    if (settings == nullptr) {
        return;
    }

    float magnitude = static_cast<float>(settings->get_setting("physics/3d/default_gravity", 9.8));
    const Vector3 direction = settings->get_setting("physics/3d/default_gravity_vector", Vector3(0, -1, 0));

    // AVBD is Z-up with gravity along -Z, so the projection of Godot's gravity vector onto Godot's
    // down direction becomes the scalar the solver wants.
    const float sign = direction.y <= 0.0f ? 1.0f : -1.0f;
    p_space.solver.gravity = -magnitude * sign;
}

double AVBDPhysicsServer3D::_space_get_param(const RID &p_space, PhysicsServer3D::SpaceParameter p_param) const {
    (void)p_space;
    (void)p_param;
    return 0.0;
}

// One state object per space, created on first ask. PhysicsDirectSpaceState3D is
// reference-counted by the engine, so the server holds a raw pointer to the object it
// owns for the space's lifetime.
PhysicsDirectSpaceState3D *AVBDPhysicsServer3D::_space_get_direct_state(const RID &p_space) {
    SpaceData *space = find_space(p_space);
    if (space == nullptr) {
        return nullptr;
    }
    if (space->direct_state == nullptr) {
        space->direct_state = memnew(AVBDDirectSpaceState3D);
        space->direct_state->bind(this, p_space);
    }
    return space->direct_state;
}

void AVBDPhysicsServer3D::_space_set_debug_contacts(const RID &p_space, int32_t p_max_contacts) {
    if (SpaceData *space = find_space(p_space)) {
        space->debug_contacts_requested = p_max_contacts > 0;
        space->debug_contacts_max = p_max_contacts;
    }
}

PackedVector3Array AVBDPhysicsServer3D::_space_get_contacts(const RID &p_space) const {
    PackedVector3Array out;
    const auto found = spaces.find(id_of(p_space));
    if (found == spaces.end()) {
        return out;
    }
    // The engine's debug contacts are world-space points of recent contacts; serve the
    // first-shape contact points of every manifold seen this step.
    for (const avbd::Force *force = found->second.solver.forces; force != nullptr; force = force->next) {
        const int count = force->contactPointCount();
        if (count <= 0) {
            continue;
        }
        // Walk the manifold's stored contacts: rA/rB are local offsets from each body.
        const avbd::Rigid *a = force->bodyA;
        const avbd::Manifold *manifold = static_cast<const avbd::Manifold *>(force);
        for (int i = 0; i < count; i++) {
            const avbd::float3 world = a->positionLin + avbd::rotate(a->positionAng, manifold->contacts[i].rA);
            out.push_back(to_godot(world));
            if (out.size() >= static_cast<int64_t>(found->second.debug_contacts_max)) {
                return out;
            }
        }
    }
    return out;
}

int32_t AVBDPhysicsServer3D::_space_get_contact_count(const RID &p_space) const {
    const auto found = spaces.find(id_of(p_space));
    if (found == spaces.end()) {
        return 0;
    }
    int32_t count = 0;
    for (const avbd::Force *force = found->second.solver.forces; force != nullptr; force = force->next) {
        count += force->contactPointCount();
    }
    return count;
}

// -----------------------------------------------------------------------------
// Shapes
//
// The engine sends a shape's parameters as a Variant from the matching *Shape3D resource, which
// is also what it can read back, so the getters rebuild the same object.
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_box_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_BOX;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_sphere_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_SPHERE;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_cylinder_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_CYLINDER;
    return make_rid(id);
}

void AVBDPhysicsServer3D::_shape_set_data(const RID &p_shape, const Variant &p_data) {
    ShapeData *shape = find_shape(p_shape);
    if (shape == nullptr) {
        return;
    }
    // Godot 4.x payload per shape type (scene/resources/3d/*_shape_3d.cpp::_update_shape):
    //   box      Vector3 - HALF extents (the resource sends `size / 2`)
    //   sphere   float radius
    //   cylinder Dictionary { radius: float, height: float } - full height
    switch (shape->type) {
        case PhysicsServer3D::SHAPE_SPHERE: {
            if (p_data.get_type() == Variant::FLOAT) {
                shape->radius = p_data;
            } else if (p_data.get_type() == Variant::DICTIONARY && ((const Dictionary)p_data).has("radius")) {
                shape->radius = ((const Dictionary)p_data)["radius"];
            }
            break;
        }
        case PhysicsServer3D::SHAPE_CYLINDER: {
            if (p_data.get_type() == Variant::DICTIONARY) {
                const Dictionary data = p_data;
                if (data.has("radius")) {
                    shape->radius = data["radius"];
                }
                if (data.has("height")) {
                    shape->height = data["height"];
                }
            }
            break;
        }
        default: {
            if (p_data.get_type() == Variant::VECTOR3) {
                // Half extents in, full extents stored (the resource's `size / 2`).
                const Vector3 half = p_data;
                shape->extents = half * 2.0;
            } else if (p_data.get_type() == Variant::DICTIONARY && ((const Dictionary)p_data).has("size")) {
                shape->extents = ((const Dictionary)p_data)["size"];
            }
            break;
        }
    }
}

PhysicsServer3D::ShapeType AVBDPhysicsServer3D::_shape_get_type(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    return found == shapes.end() ? PhysicsServer3D::SHAPE_BOX : found->second.type;
}

Variant AVBDPhysicsServer3D::_shape_get_data(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    if (found == shapes.end()) {
        return Variant();
    }
    const ShapeData &shape = found->second;

    // The exact form the engine's resources send: box = half extents Vector3,
    // sphere = radius float, cylinder = Dictionary.
    switch (shape.type) {
        case PhysicsServer3D::SHAPE_SPHERE:
            return Variant(shape.radius);
        case PhysicsServer3D::SHAPE_CYLINDER: {
            Dictionary data;
            data["radius"] = shape.radius;
            data["height"] = shape.height;
            return data;
        }
        default:
            return Variant(shape.extents / 2.0);
    }
}

// Godot BodyAxis bits (linear X/Y/Z = 1/2/4, angular X/Y/Z = 8/16/32) to the solver's
// sim-axis bits (0..2 = sim x/y/z). Godot Y is sim Z, Godot Z is sim -Y; the lock is a
// directionless bit so the sign does not matter, only the permutation.
static void axis_locks_to_sim(uint32_t p_locks, uint8_t &r_linear, uint8_t &r_angular) {
    static constexpr int kMap[3] = {0, 2, 1}; // godot axis -> sim axis
    r_linear = 0;
    r_angular = 0;
    for (int godot_axis = 0; godot_axis < 3; godot_axis++) {
        if (p_locks & (1u << godot_axis)) {
            r_linear |= static_cast<uint8_t>(1u << kMap[godot_axis]);
        }
        if (p_locks & (1u << (godot_axis + 3))) {
            r_angular |= static_cast<uint8_t>(1u << kMap[godot_axis]);
        }
    }
}

// -----------------------------------------------------------------------------
// The solver
//
// A space's solver is rebuilt from its bodies whenever any of them changes: AVBD holds one shape
// per body and a flat body list, so there is nothing to incrementally patch. Bodies are few
// enough per space that rebuilding is cheap, and it keeps the solver a pure function of the scene
// rather than of the order the engine happened to call things in.
// -----------------------------------------------------------------------------

// Density the solver should see for a body: its mass divided by its first shape's volume. The
// solver derives mass and inertia from shape and density, so the density has to carry exactly
// the engine's mass - a body that never set one is 1 kg (BODY_PARAM_MASS's documented default),
// whatever its size. Mass 0 means a body with no volume to divide into, or one the engine
// explicitly made massless; both stay massless in the solver, which reads that as immovable.
float AVBDPhysicsServer3D::density_for(const BodyData &p_body) {
    if (p_body.mode != PhysicsServer3D::BODY_MODE_RIGID &&
            p_body.mode != PhysicsServer3D::BODY_MODE_RIGID_LINEAR) {
        return 0.0f; // static and kinematic bodies are immovable
    }
    if (p_body.shapes.empty() || p_body.mass <= 0.0) {
        return 0.0f;
    }
    const ShapeData &shape = p_body.shapes.front();
    double volume = 0.0;
    switch (shape.type) {
        case PhysicsServer3D::SHAPE_SPHERE:
            volume = (4.0 / 3.0) * 3.14159265358979 * shape.radius * shape.radius * shape.radius;
            break;
        case PhysicsServer3D::SHAPE_CYLINDER:
            volume = 3.14159265358979 * shape.radius * shape.radius * shape.height;
            break;
        default:
            volume = shape.extents.x * shape.extents.y * shape.extents.z;
            break;
    }
    if (volume <= 1.0e-9) {
        return 0.0f;
    }
    return static_cast<float>(p_body.mass / volume);
}

void AVBDPhysicsServer3D::rebuild_space(const RID &p_space) {
    SpaceData *space = find_space(p_space);
    if (space == nullptr) {
        return;
    }

    const uint64_t space_id = id_of(p_space);

    // Anything this space owned is gone; the solver is about to be refilled.
    for (auto &[id, body] : bodies) {
        if (id_of(body.space) == space_id) {
            body.rigid = nullptr;
        }
    }
    space->solver.clear();

    // Deterministic order: ids are minted in creation order, so sorting by id makes the
    // solver's body list - and with it the Gauss-Seidel update order - a function of the
    // scene, not of unordered_map bucket layout. Two identical scenes rebuild identically.
    std::vector<uint64_t> space_body_ids;
    for (const auto &[id, body] : bodies) {
        if (id_of(body.space) == space_id) {
            space_body_ids.push_back(id);
        }
    }
    std::sort(space_body_ids.begin(), space_body_ids.end());

    for (const uint64_t id : space_body_ids) {
        BodyData &body = bodies.find(id)->second;

        // Refresh the shape snapshots from the engine's shape records: the resource data
        // may have arrived after the shape was attached to the body.
        for (size_t i = 0; i < body.shapes.size() && i < body.shape_rids.size(); i++) {
            if (const ShapeData *fresh = shape_data(body.shape_rids[i])) {
                body.shapes[i] = *fresh;
            }
        }

        // A body with no shape cannot collide and has no inertia; the engine always attaches one
        // before it matters, and until then the body is left out of the solver.
        if (body.shapes.empty()) {
            continue;
        }

        const ShapeData &shape = body.shapes.front();
        if (body.shapes.size() > 1) {
            WARN_PRINT("AVBD physics server: a body has more than one shape; AVBD simulates one "
                       "shape per body, so the first is used.");
        }

        avbd::float3 extents{1, 1, 1};
        avbd::ShapeType type = avbd::ShapeType::Box;
        switch (shape.type) {
            case PhysicsServer3D::SHAPE_SPHERE:
                type = avbd::ShapeType::Sphere;
                extents = {static_cast<float>(shape.radius), 0, 0};
                break;
            case PhysicsServer3D::SHAPE_CYLINDER:
                type = avbd::ShapeType::Cylinder;
                // Godot: x and z are the radius, y the height, axis local Y. The solver reads
                // the radius from size.x and the height from size.z (its cylinder axis is its
                // local Z, which is where Godot's Y lands in sim space), so order {r, r, h}.
                extents = {static_cast<float>(shape.radius), static_cast<float>(shape.radius),
                        static_cast<float>(shape.height)};
                break;
            default:
                // Godot's BoxShape3D size is full extents, and so is the solver's box size;
                // permute into sim axes rather than scaling.
                extents = to_sim_extents(shape.extents);
                break;
        }

        body.rigid = new avbd::Rigid(&space->solver, extents, type, density_for(body), body.friction,
                to_sim(body.transform.origin), to_sim(body.linear_velocity));
        // The same Y-up to Z-up mapping the node layer uses, applied to the orientation too.
        body.rigid->positionAng = to_sim(body.transform.basis.get_rotation_quaternion());
        body.rigid->velocityAng = to_sim(body.angular_velocity);
        body.rigid->collisionLayer = body.collision_layer;
        body.rigid->collisionMask = body.collision_mask;
        axis_locks_to_sim(body.axis_locks, body.rigid->axisLockLinear, body.rigid->axisLockAngular);
        body.rigid->sleeping = body.sleeping;
        // Solver sleep modes: 0 = never sleeps, 1 = may sleep. Sleep is manual this round:
        // the solver has no idle detection, so "can sleep" only opens the gate.
        body.rigid->sleep_mode = body.can_sleep ? 1 : 0;
        body.rigid->gravity = space->solver.gravity * static_cast<float>(body.gravity_scale);
    }

    // Collision exceptions: a pair never collides while either side lists the other. One
    // solver IgnoreCollision per pair; its mere existence is what the broad phase tests.
    for (auto &[id, body] : bodies) {
        if (id_of(body.space) != space_id || body.rigid == nullptr) {
            continue;
        }
        for (const uint64_t other_id : body.exceptions) {
            if (other_id <= id) {
                continue; // each unordered pair once (from its lower id)
            }
            const auto other_found = bodies.find(other_id);
            if (other_found == bodies.end() || other_found->second.rigid == nullptr ||
                    id_of(other_found->second.space) != space_id) {
                continue;
            }
            new avbd::IgnoreCollision(&space->solver, body.rigid, other_found->second.rigid);
        }
    }

    // Joints second, once every body exists. One solver GenericJoint per engine joint
    // whose both bodies live in this space (or whose A side is the world). Every Godot
    // joint type maps onto the solver's GenericJoint; the specialised types fix which
    // axes are locked (conventions verified against Godot 4.7 sources):
    //   Pin       - all three linear axes locked, all angular free
    //   Hinge     - linear locked; angular free about the hinge axis, locked elsewhere
    //               (the hinge axis is the Z column of the hinge frame, see
    //               modules/godot_physics_3d/joints/godot_hinge_joint_3d.cpp:160)
    //   Slider    - free (or limited) along the slide axis (the X column of the joint
    //               frame, godot_slider_joint_3d.cpp:125), locked on the orthogonal
    //               linear axes; angular free about the slide axis, locked elsewhere
    //   ConeTwist - linear locked; angular free about the twist axis, swing limited
    //               on the two orthogonal axes (godot_cone_twist_joint_3d.cpp:117-139)
    // A joint axis that does not line up with a solver axis snaps to the nearest one:
    // the GenericJoint's degrees of freedom are principal axes only.
    std::vector<uint64_t> space_joint_ids;
    for (const auto &[id, joint] : joints) {
        (void)joint;
        space_joint_ids.push_back(id);
    }
    std::sort(space_joint_ids.begin(), space_joint_ids.end());
    for (const uint64_t joint_id : space_joint_ids) {
        const JointData &joint = joints.find(joint_id)->second;
        const BodyData *data_a = joint.body_a.is_valid() ? find_body(joint.body_a) : nullptr;
        const BodyData *data_b = find_body(joint.body_b);
        if (data_b == nullptr || data_b->rigid == nullptr) {
            continue;
        }
        if (joint.body_a.is_valid() && (data_a == nullptr || data_a->rigid == nullptr)) {
            continue;
        }
        avbd::Rigid *rigid_a = joint.body_a.is_valid() ? data_a->rigid : nullptr;
        if (rigid_a == nullptr && id_of(data_b->space) != space_id) {
            continue;
        }
        if (rigid_a != nullptr && id_of(data_a->space) != space_id) {
            continue;
        }
        if (joint.type == PhysicsServer3D::JOINT_TYPE_6DOF) {
            rebuild_6dof_joint(space, joint, rigid_a, data_b->rigid);
            continue;
        }
        if (joint.type == PhysicsServer3D::JOINT_TYPE_PIN) {
            rebuild_pin_joint(space, joint, rigid_a, data_b->rigid);
            continue;
        }
        if (joint.type == PhysicsServer3D::JOINT_TYPE_HINGE) {
            rebuild_hinge_joint(space, joint, rigid_a, data_b->rigid);
            continue;
        }
        if (joint.type == PhysicsServer3D::JOINT_TYPE_SLIDER) {
            rebuild_slider_joint(space, joint, rigid_a, data_b->rigid);
            continue;
        }
        if (joint.type == PhysicsServer3D::JOINT_TYPE_CONE_TWIST) {
            rebuild_cone_twist_joint(space, joint, rigid_a, data_b->rigid);
            continue;
        }
    }

    // Constant forces that were queued before this body existed in the solver are applied on
    // the first step (see _step); nothing to do here.
}

// --- specialised joint mappings -------------------------------------------------

// The solver sim axis best aligned with a (sim-space) joint axis. The GenericJoint's
// degrees of freedom are principal axes, so an off-axis hinge snaps to the dominant one.
static int dominant_sim_axis(const avbd::float3 &axis_sim) {
    const float abs_axis[3] = {std::fabs(axis_sim.x), std::fabs(axis_sim.y), std::fabs(axis_sim.z)};
    int best = 0;
    if (abs_axis[1] > abs_axis[best]) {
        best = 1;
    }
    if (abs_axis[2] > abs_axis[best]) {
        best = 2;
    }
    return best;
}

// Sign of the axis' dominant component (+1 / -1), for mirroring limits.
static float dominant_sim_axis_sign(const avbd::float3 &axis_sim, int solver_axis) {
    const float value = solver_axis == 0 ? axis_sim.x : (solver_axis == 1 ? axis_sim.y : axis_sim.z);
    return value >= 0.0f ? 1.0f : -1.0f;
}

// Const-safe lookups into a joint record's parameter maps (default 0 when unset).
static double param_value(const AVBDPhysicsServer3D::JointData &p_joint, int p_key) {
    const auto found = p_joint.joint_params.find(p_key);
    return found == p_joint.joint_params.end() ? 0.0 : found->second;
}

static bool param_flag(const AVBDPhysicsServer3D::JointData &p_joint, int p_key) {
    const auto found = p_joint.joint_flags.find(p_key);
    return found != p_joint.joint_flags.end() && found->second;
}

void AVBDPhysicsServer3D::rebuild_6dof_joint(SpaceData *space, const JointData &joint, avbd::Rigid *rigid_a,
        avbd::Rigid *rigid_b) {
    // Anchors: Godot's local_ref is a full transform whose basis carries the joint frame;
    // AVBD wants the offset in each body's frame and measures angles from the relative
    // orientation at creation, which is what the two bases encode. With no A body the A
    // anchor is world-space already (the server's sim space is the world).
    const avbd::float3 r_a = to_sim(joint.ref_a.origin);
    const avbd::float3 r_b = to_sim(joint.ref_b.origin);

    avbd::GenericJoint *solver_joint = new avbd::GenericJoint(&space->solver, rigid_a,
            rigid_b, r_a, r_b);

    // Six axes, Godot order X..Z linear then angular. Each Godot axis maps to a solver axis
    // with a sign (Godot X -> solver X +, Y -> Z +, Z -> Y -); the same kAxisMapping the
    // node layer uses, so both frontends mean the same joint.
    for (int axis = 0; axis < 6; axis++) {
        const int godot_index = axis % 3;
        const bool angular = axis >= 3;
        const AxisMapping &map = kAxisMapping[godot_index];
        avbd::JointAxis &target = angular ? solver_joint->angular[map.solver]
                                          : solver_joint->linear[map.solver];

        const double lower = joint.params[axis][PhysicsServer3D::G6DOF_JOINT_LINEAR_LOWER_LIMIT];
        const double upper = joint.params[axis][PhysicsServer3D::G6DOF_JOINT_LINEAR_UPPER_LIMIT];
        const bool limit_enabled = angular
                ? joint.flags[axis][PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_ANGULAR_LIMIT]
                : joint.flags[axis][PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_LINEAR_LIMIT];

        if (!limit_enabled) {
            target.mode = avbd::AxisMode::Free;
        } else if (std::fabs(upper - lower) < 1.0e-9) {
            target.mode = avbd::AxisMode::Locked; // Godot spells "locked" as limits 0..0
        } else {
            target.mode = avbd::AxisMode::Limited;
        }

        // Springs: stiffness above zero and the spring flag on (Godot 4.3 enum slots;
        // values at the top of this file since the pinned godot-cpp headers lack them).
        const int spring_flag = angular ? g6dof_spring::FLAG_ENABLE_ANGULAR
                                        : g6dof_spring::FLAG_ENABLE_LINEAR;
        const double stiffness = angular ? joint.params[axis][g6dof_spring::ANGULAR_STIFFNESS]
                                         : joint.params[axis][g6dof_spring::LINEAR_STIFFNESS];
        const double damping = angular ? joint.params[axis][g6dof_spring::ANGULAR_DAMPING]
                                       : joint.params[axis][g6dof_spring::LINEAR_DAMPING];
        const double equilibrium = angular
                ? joint.params[axis][g6dof_spring::ANGULAR_EQUILIBRIUM]
                : joint.params[axis][g6dof_spring::LINEAR_EQUILIBRIUM];
        if (joint.flags[axis][spring_flag] && stiffness > 0.0) {
            target.springStiffness = static_cast<float>(stiffness);
            target.springDamping = static_cast<float>(damping);
            target.springEquilibrium = static_cast<float>(equilibrium) * map.sign;
        }

        // A reversed axis mirrors its limits and equilibrium, as in the node layer.
        if (map.sign > 0.0f) {
            target.lower = static_cast<float>(lower);
            target.upper = static_cast<float>(upper);
        } else {
            target.lower = static_cast<float>(-upper);
            target.upper = static_cast<float>(-lower);
        }
    }
}

void AVBDPhysicsServer3D::rebuild_pin_joint(SpaceData *space, const JointData &joint, avbd::Rigid *rigid_a,
        avbd::Rigid *rigid_b) {
    const avbd::float3 r_a = to_sim(rigid_a != nullptr ? joint.pin_local_a : joint.pin_local_a);
    const avbd::float3 r_b = to_sim(joint.pin_local_b);

    avbd::GenericJoint *solver_joint = new avbd::GenericJoint(&space->solver, rigid_a, rigid_b, r_a, r_b);

    // A pin is a ball joint: all three linear axes locked, rotation free.
    for (int solver_axis = 0; solver_axis < 3; solver_axis++) {
        solver_joint->linear[solver_axis].mode = avbd::AxisMode::Locked;
    }
    for (int solver_axis = 0; solver_axis < 3; solver_axis++) {
        solver_joint->angular[solver_axis].mode = avbd::AxisMode::Free;
    }
    // PIN_JOINT_BIAS / DAMPING / IMPULSE_CLAMP: recorded, not mapped (the solver's
    // stabilisation is its own; documented in the README).
}

void AVBDPhysicsServer3D::rebuild_hinge_joint(SpaceData *space, const JointData &joint, avbd::Rigid *rigid_a,
        avbd::Rigid *rigid_b) {
    // The hinge axis: for _joint_make_hinge_simple it is the vector the engine stored in
    // ref_a's basis Z column (see _joint_make_hinge_simple below); for the full form it is
    // the Z column of the hinge frame. Expressed in body A's local frame.
    const avbd::float3 axis_a_sim = to_sim(joint.ref_a.basis.get_column(2));
    const int hinge_axis = dominant_sim_axis(axis_a_sim);
    const float hinge_sign = dominant_sim_axis_sign(axis_a_sim, hinge_axis);

    const avbd::float3 r_a = to_sim(joint.ref_a.origin);
    const avbd::float3 r_b = to_sim(joint.ref_b.origin);

    avbd::GenericJoint *solver_joint = new avbd::GenericJoint(&space->solver, rigid_a, rigid_b, r_a, r_b);

    // Godot's hinge limits sit around the joint's rest angle with lower <= 0 <= upper.
    const bool use_limit = param_flag(joint, PhysicsServer3D::HINGE_JOINT_FLAG_USE_LIMIT);
    const double upper = param_value(joint, PhysicsServer3D::HINGE_JOINT_LIMIT_UPPER);
    const double lower = param_value(joint, PhysicsServer3D::HINGE_JOINT_LIMIT_LOWER);

    for (int solver_axis = 0; solver_axis < 3; solver_axis++) {
        solver_joint->linear[solver_axis].mode = avbd::AxisMode::Locked;
        if (solver_axis == hinge_axis) {
            if (use_limit && std::fabs(upper - lower) > 1.0e-9) {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Limited;
                solver_joint->angular[solver_axis].lower = static_cast<float>(hinge_sign > 0.0f ? lower : -upper);
                solver_joint->angular[solver_axis].upper = static_cast<float>(hinge_sign > 0.0f ? upper : -lower);
            } else {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Free;
            }
        } else {
            solver_joint->angular[solver_axis].mode = avbd::AxisMode::Locked;
        }
    }
    // HINGE_JOINT_FLAG_ENABLE_MOTOR and its velocity/force params: recorded, not mapped.
}

void AVBDPhysicsServer3D::rebuild_slider_joint(SpaceData *space, const JointData &joint, avbd::Rigid *rigid_a,
        avbd::Rigid *rigid_b) {
    // The slide axis is the X column of the joint frame (godot_slider_joint_3d.cpp:125).
    const avbd::float3 axis_a_sim = to_sim(joint.ref_a.basis.get_column(0));
    const int slide_axis = dominant_sim_axis(axis_a_sim);
    const float slide_sign = dominant_sim_axis_sign(axis_a_sim, slide_axis);

    const avbd::float3 r_a = to_sim(joint.ref_a.origin);
    const avbd::float3 r_b = to_sim(joint.ref_b.origin);

    avbd::GenericJoint *solver_joint = new avbd::GenericJoint(&space->solver, rigid_a, rigid_b, r_a, r_b);

    const double lin_upper = param_value(joint, PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_UPPER);
    const double lin_lower = param_value(joint, PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_LOWER);
    const double ang_upper = param_value(joint, PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_UPPER);
    const double ang_lower = param_value(joint, PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_LOWER);

    for (int solver_axis = 0; solver_axis < 3; solver_axis++) {
        // Linear: slide axis free or limited, the orthogonal pair locked.
        if (solver_axis == slide_axis) {
            if (std::fabs(lin_upper - lin_lower) > 1.0e-9) {
                solver_joint->linear[solver_axis].mode = avbd::AxisMode::Limited;
                solver_joint->linear[solver_axis].lower = static_cast<float>(slide_sign > 0.0f ? lin_lower : -lin_upper);
                solver_joint->linear[solver_axis].upper = static_cast<float>(slide_sign > 0.0f ? lin_upper : -lin_lower);
            } else {
                solver_joint->linear[solver_axis].mode = avbd::AxisMode::Locked;
            }
        } else {
            solver_joint->linear[solver_axis].mode = avbd::AxisMode::Locked;
        }

        // Angular: rotation about the slide axis limited per the angular limits
        // (defaults 0/0 = locked in Godot), the orthogonal pair locked.
        if (solver_axis == slide_axis) {
            if (std::fabs(ang_upper - ang_lower) > 1.0e-9) {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Limited;
                solver_joint->angular[solver_axis].lower = static_cast<float>(slide_sign > 0.0f ? ang_lower : -ang_upper);
                solver_joint->angular[solver_axis].upper = static_cast<float>(slide_sign > 0.0f ? ang_upper : -ang_lower);
            } else {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Locked;
            }
        } else {
            solver_joint->angular[solver_axis].mode = avbd::AxisMode::Locked;
        }
    }
    // The softness / restitution / damping / motor params of the slider family: recorded,
    // not mapped (documented).
}

void AVBDPhysicsServer3D::rebuild_cone_twist_joint(SpaceData *space, const JointData &joint, avbd::Rigid *rigid_a,
        avbd::Rigid *rigid_b) {
    // The twist axis is the X column of the joint frame (b1Axis1 at
    // godot_cone_twist_joint_3d.cpp:117); the two orthogonal axes carry the swing span.
    const avbd::float3 axis_a_sim = to_sim(joint.ref_a.basis.get_column(0));
    const int twist_axis = dominant_sim_axis(axis_a_sim);
    const float twist_sign = dominant_sim_axis_sign(axis_a_sim, twist_axis);

    const avbd::float3 r_a = to_sim(joint.ref_a.origin);
    const avbd::float3 r_b = to_sim(joint.ref_b.origin);

    avbd::GenericJoint *solver_joint = new avbd::GenericJoint(&space->solver, rigid_a, rigid_b, r_a, r_b);

    const double swing = param_value(joint, PhysicsServer3D::CONE_TWIST_JOINT_SWING_SPAN);
    const double twist = param_value(joint, PhysicsServer3D::CONE_TWIST_JOINT_TWIST_SPAN);

    for (int solver_axis = 0; solver_axis < 3; solver_axis++) {
        solver_joint->linear[solver_axis].mode = avbd::AxisMode::Locked;
        if (solver_axis == twist_axis) {
            // Twist span limits rotation about the twist axis (Godot's default is 0,
            // i.e. no twist; a negative span means unlimited).
            if (twist < -1.0e-9) {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Free;
            } else if (twist > 1.0e-9) {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Limited;
                solver_joint->angular[solver_axis].lower = static_cast<float>(twist_sign > 0.0f ? -twist : -twist);
                solver_joint->angular[solver_axis].upper = static_cast<float>(twist);
            } else {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Locked;
            }
        } else {
            // Swing: Godot applies the same span to both orthogonal axes, centred.
            if (swing > 1.0e-9) {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Limited;
                solver_joint->angular[solver_axis].lower = static_cast<float>(-swing * 0.5f);
                solver_joint->angular[solver_axis].upper = static_cast<float>(swing * 0.5f);
            } else {
                solver_joint->angular[solver_axis].mode = avbd::AxisMode::Locked;
            }
        }
    }
    // CONE_TWIST_JOINT_BIAS / SOFTNESS / RELAXATION: recorded, not mapped.
}

// Run every active space for one step, then tell each body what happened.
//
// Stepping follows the game: while the engine has the server inactive (SceneTree paused)
// nothing steps. Solver parameters are stamped from project settings every step, and each
// space runs `substeps` solver steps per engine step with dt split across them.
void AVBDPhysicsServer3D::_step(double p_step) {
    if (!server_active) {
        return;
    }

    // One state object serves every callback this step; re-aimed per body before each call.
    if (body_state == nullptr) {
        body_state = memnew(AVBDDirectBodyState3D);
    }

    for (auto &[id, space] : spaces) {
        if (!space.active) {
            continue;
        }

        // A body joined/left this space since the last step; rebuild now that the
        // engine has finished sending its setup (shapes, mass, material). Runs before
        // the empty check: a space may be pending precisely because its first body is
        // only now arriving.
        if (space_rebuild_pending.erase(id) > 0) {
            rebuild_space(make_rid(id));
        }

        if (space.solver.bodies == nullptr) {
            continue;
        }

        // Project settings may have changed since the last step (the tests flip
        // physics/avbd/substeps at runtime); re-read before stamping.
        read_project_params(space);

        const float dt = static_cast<float>(p_step);
        const int substeps = space.substeps > 0 ? space.substeps : 1;
        space.solver.dt = dt / static_cast<float>(substeps);
        space.solver.iterations = space.iterations > 0 ? space.iterations : space.solver.iterations;
        space.solver.alpha = static_cast<float>(space.alpha);
        space.solver.betaLin = static_cast<float>(space.beta_linear);
        space.solver.betaAng = static_cast<float>(space.beta_angular);
        space.solver.gamma = static_cast<float>(space.gamma);
        space.solver.threads = space.threads;

        // Forces accumulated since the last step (persistent constant_* plus one-shot frame
        // ones) become velocity changes now, so the solver integrates them with everything else.
        for (auto &[body_id, body] : bodies) {
            if (id_of(body.space) != id || body.rigid == nullptr || body.rigid->mass <= 0.0f) {
                continue;
            }
            // A body with scripted integration applies its own forces in its callback.
            if (body.omit_force_integration) {
                body.frame_force = Vector3();
                body.frame_torque = Vector3();
                continue;
            }
            const Vector3 force = body.constant_force + body.frame_force;
            const Vector3 torque = body.constant_torque + body.frame_torque;
            body.frame_force = Vector3();
            body.frame_torque = Vector3();
            if (force.length_squared() > 0.0) {
                body.rigid->velocityLin += to_sim(force) / body.rigid->mass * dt;
            }
            if (torque.length_squared() > 0.0) {
                // World-frame inverse inertia, as in the node layer's impulse paths.
                const avbd::float3 local = avbd::rotate(avbd::conjugate(body.rigid->positionAng),
                        to_sim(torque) * dt);
                const avbd::float3 delta{local.x / body.rigid->moment.x,
                        local.y / body.rigid->moment.y, local.z / body.rigid->moment.z};
                body.rigid->velocityAng += avbd::rotate(body.rigid->positionAng, delta);
            }
        }

        for (int sub = 0; sub < substeps; sub++) {
            space.solver.step();
        }

        // Area monitoring runs once per engine step, after the solver's output is in hand
        // and before the scene reads anything back.
        _monitor_areas();

        // Collect each body's contacts for the step (used by the DirectBodyState contact
        // getters, i.e. RigidBody3D contact monitoring). Manifolds are the forces with
        // contact points; a contact's world position comes from the anchor of the body
        // being filled. Capped by max_contacts_reported, as Godot does: without the cap
        // request no contacts are recorded.
        for (auto &[collect_id, collect_body] : bodies) {
            (void)collect_id;
            if (id_of(collect_body.space) != id) {
                collect_body.frame_contacts.clear();
            }
        }
        {
            // Rigid -> owning record, resolved once for the whole manifold walk.
            std::unordered_map<const avbd::Rigid *, BodyData *> rigid_map;
            std::unordered_map<const BodyData *, uint64_t> id_map;
            for (auto &[bid, b] : bodies) {
                if (id_of(b.space) == id && b.rigid != nullptr) {
                    b.frame_contacts.clear();
                    rigid_map.emplace(b.rigid, &b);
                    id_map.emplace(&b, bid);
                }
            }
            for (avbd::Force *force = space.solver.forces; force != nullptr; force = force->next) {
                const int count = force->contactPointCount();
                if (count <= 0) {
                    continue;
                }
                const avbd::Manifold *manifold = static_cast<const avbd::Manifold *>(force);
                const auto it_a = rigid_map.find(force->bodyA);
                const auto it_b = rigid_map.find(force->bodyB);
                BodyData *data[2] = {it_a != rigid_map.end() ? it_a->second : nullptr,
                        it_b != rigid_map.end() ? it_b->second : nullptr};

                for (int side = 0; side < 2; side++) {
                    BodyData *fill = data[side];
                    if (fill == nullptr || fill->max_contacts_reported <= 0) {
                        continue;
                    }
                    const BodyData *other = data[side == 0 ? 1 : 0];
                    const avbd::Rigid *self_rigid = side == 0 ? force->bodyA : force->bodyB;
                    for (int c = 0; c < count && static_cast<int32_t>(fill->frame_contacts.size()) < fill->max_contacts_reported; c++) {
                        const avbd::float3 world = self_rigid->positionLin
                                + avbd::rotate(self_rigid->positionAng,
                                        side == 0 ? manifold->contacts[c].rA : manifold->contacts[c].rB);
                        BodyData::FrameContact fc;
                        fc.position = to_godot(world);
                        // The basis's first row points from B to A; each side reports the
                        // normal as pointing out of the other body towards itself.
                        fc.normal = to_godot(manifold->basis[0] * (side == 0 ? 1.0f : -1.0f));
                        fc.collider_body_id = other != nullptr ? id_map[other] : 0;
                        fill->frame_contacts.push_back(fc);
                    }
                }
            }
        }

        for (auto &[body_id, body] : bodies) {
            if (id_of(body.space) != id || body.rigid == nullptr) {
                continue;
            }

            // Wall clock the result back, through the callback the engine registered for this
            // body. The server does not touch the scene itself: the engine applies it.
            const avbd::Rigid *rigid = body.rigid;

            Transform3D transform;
            transform.origin = to_godot(rigid->positionLin);
            transform.basis = Basis(to_godot(rigid->positionAng).normalized());

            body.transform = transform;
            body.linear_velocity = to_godot(rigid->velocityLin);
            body.angular_velocity = to_godot(rigid->velocityAng);

            // Godot 4.x sync contract (verified against 4.7 sources): the callable receives
            // exactly one argument - the PhysicsDirectBodyState3D. The engine reads transform,
            // velocities, sleeping and contacts off it.
            if (body.sync.is_valid()) {
                body_state->bind(this, body.rigid, make_rid(body_id));
                Array args;
                args.push_back(body_state);
                body.sync.callv(args);
            }

            // Scripted integration: the callback gets the state plus the stored userdata
            // (Godot drops the userdata argument when it is nil).
            if (body.force_integration_callback.is_valid()) {
                body_state->bind(this, body.rigid, make_rid(body_id));
                Array args;
                args.push_back(body_state);
                const bool has_userdata = body.force_integration_userdata.get_type() != Variant::NIL;
                if (has_userdata) {
                    args.push_back(body.force_integration_userdata);
                }
                body.force_integration_callback.callv(args);
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Joints
//
// The engine hands over a joint RID, then a type (joint_make_*), then axis parameters and flags
// one call at a time. Everything is recorded in JointData; the solver object is built only in
// rebuild_space, so a joint whose bodies are still being set up costs nothing.
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_joint_create() {
    const uint64_t id = next_id++;
    joints.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple());
    return make_rid(id);
}

void AVBDPhysicsServer3D::_joint_clear(const RID &p_joint) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    // Clearing resets the record to an untyped joint, as the engine's contract states; bodies
    // the joint constrained get a rebuild so the solver constraint goes away.
    const RID body_a = found->second.body_a;
    const RID body_b = found->second.body_b;
    found->second = JointData();
    rebuild_space(find_body(body_a) ? find_body(body_a)->space : RID());
    if (BodyData *b = find_body(body_b)) {
        rebuild_space(b->space);
    }
}

PhysicsServer3D::JointType AVBDPhysicsServer3D::_joint_get_type(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? PhysicsServer3D::JOINT_TYPE_6DOF : found->second.type;
}

void AVBDPhysicsServer3D::_joint_set_solver_priority(const RID &p_joint, int32_t p_priority) {
    const auto found = joints.find(id_of(p_joint));
    if (found != joints.end()) {
        found->second.priority = p_priority; // AVBD has one solver; ordering is not modelled
    }
}

int32_t AVBDPhysicsServer3D::_joint_get_solver_priority(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? 1 : found->second.priority;
}

void AVBDPhysicsServer3D::_joint_disable_collisions_between_bodies(const RID &p_joint, bool p_disable) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    // A solver joint keeps its bodies from colliding no matter what (they share a force), so
    // Godot's default exclusion comes free; asking for contact between the pair cannot be
    // honoured and the record only reflects the request.
    found->second.exclude_collision = p_disable;
}

bool AVBDPhysicsServer3D::_joint_is_disabled_collisions_between_bodies(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? true : found->second.exclude_collision;
}

void AVBDPhysicsServer3D::_joint_make_generic_6dof(const RID &p_joint, const RID &p_body_a,
        const Transform3D &p_local_ref_a, const RID &p_body_b, const Transform3D &p_local_ref_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_6DOF;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_local_ref_a;
    joint.ref_b = p_local_ref_b;

    if (BodyData *body_b = find_body(p_body_b)) {
        rebuild_space(body_b->space);
    }
}

// A shared body lookup that triggers the rebuild a completed joint definition needs.
void AVBDPhysicsServer3D::_rebuild_for_joint(const RID &p_joint) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    if (BodyData *body_b = find_body(found->second.body_b)) {
        rebuild_space(body_b->space);
    } else if (BodyData *body_a = find_body(found->second.body_a)) {
        rebuild_space(body_a->space);
    }
}

void AVBDPhysicsServer3D::_joint_make_pin(const RID &p_joint, const RID &p_body_a, const Vector3 &p_local_a,
        const RID &p_body_b, const Vector3 &p_local_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_PIN;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.pin_local_a = p_local_a;
    joint.pin_local_b = p_local_b;
    joint.ref_a = Transform3D(Basis(), p_local_a);
    joint.ref_b = Transform3D(Basis(), p_local_b);
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_pin_joint_set_param(const RID &p_joint, PhysicsServer3D::PinJointParam p_param,
        double p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    // BIAS/DAMPING/IMPULSE_CLAMP are recorded only; the solver's stabilisation is its own.
}

double AVBDPhysicsServer3D::_pin_joint_get_param(const RID &p_joint, PhysicsServer3D::PinJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_pin_joint_set_local_a(const RID &p_joint, const Vector3 &p_local_a) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.pin_local_a = p_local_a;
    found->second.ref_a = Transform3D(Basis(), p_local_a);
    _rebuild_for_joint(p_joint);
}

Vector3 AVBDPhysicsServer3D::_pin_joint_get_local_a(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? Vector3() : found->second.pin_local_a;
}

void AVBDPhysicsServer3D::_pin_joint_set_local_b(const RID &p_joint, const Vector3 &p_local_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.pin_local_b = p_local_b;
    found->second.ref_b = Transform3D(Basis(), p_local_b);
    _rebuild_for_joint(p_joint);
}

Vector3 AVBDPhysicsServer3D::_pin_joint_get_local_b(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? Vector3() : found->second.pin_local_b;
}

void AVBDPhysicsServer3D::_joint_make_hinge(const RID &p_joint, const RID &p_body_a, const Transform3D &p_hinge_a,
        const RID &p_body_b, const Transform3D &p_hinge_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_HINGE;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_hinge_a;
    joint.ref_b = p_hinge_b;
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_joint_make_hinge_simple(const RID &p_joint, const RID &p_body_a,
        const Vector3 &p_pivot_a, const Vector3 &p_axis_a, const RID &p_body_b, const Vector3 &p_pivot_b,
        const Vector3 &p_axis_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_HINGE;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    // Store the hinge axis in the Z column of the stored frames, matching the full form's
    // convention (the hinge axis is the frame's Z, per godot_hinge_joint_3d.cpp:160).
    const Vector3 axis_a = p_axis_a.normalized();
    const Vector3 axis_b = p_axis_b.normalized();
    Vector3 ortho_a = std::fabs(axis_a.z) < 0.9 ? Vector3(0, 0, 1) : Vector3(1, 0, 0);
    Vector3 ortho_b = std::fabs(axis_b.z) < 0.9 ? Vector3(0, 0, 1) : Vector3(1, 0, 0);
    joint.ref_a = Transform3D(Basis(ortho_a.cross(axis_a).normalized(), axis_a.cross(ortho_a.cross(axis_a)).normalized(),
            axis_a), p_pivot_a);
    joint.ref_b = Transform3D(Basis(ortho_b.cross(axis_b).normalized(), axis_b.cross(ortho_b.cross(axis_b)).normalized(),
            axis_b), p_pivot_b);
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_hinge_joint_set_param(const RID &p_joint, PhysicsServer3D::HingeJointParam p_param,
        double p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    if (p_param == PhysicsServer3D::HINGE_JOINT_LIMIT_LOWER || p_param == PhysicsServer3D::HINGE_JOINT_LIMIT_UPPER) {
        _rebuild_for_joint(p_joint);
    }
}

double AVBDPhysicsServer3D::_hinge_joint_get_param(const RID &p_joint, PhysicsServer3D::HingeJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_hinge_joint_set_flag(const RID &p_joint, PhysicsServer3D::HingeJointFlag p_flag,
        bool p_enabled) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_flags[static_cast<int>(p_flag)] = p_enabled;
    if (p_flag == PhysicsServer3D::HINGE_JOINT_FLAG_USE_LIMIT) {
        _rebuild_for_joint(p_joint);
    }
}

bool AVBDPhysicsServer3D::_hinge_joint_get_flag(const RID &p_joint, PhysicsServer3D::HingeJointFlag p_flag) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return false;
    }
    const auto flag = found->second.joint_flags.find(static_cast<int>(p_flag));
    return flag != found->second.joint_flags.end() && flag->second;
}

void AVBDPhysicsServer3D::_joint_make_slider(const RID &p_joint, const RID &p_body_a,
        const Transform3D &p_local_ref_a, const RID &p_body_b, const Transform3D &p_local_ref_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_SLIDER;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_local_ref_a;
    joint.ref_b = p_local_ref_b;
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_slider_joint_set_param(const RID &p_joint, PhysicsServer3D::SliderJointParam p_param,
        double p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    if (p_param == PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_UPPER ||
            p_param == PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_LOWER ||
            p_param == PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_UPPER ||
            p_param == PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_LOWER) {
        _rebuild_for_joint(p_joint);
    }
}

double AVBDPhysicsServer3D::_slider_joint_get_param(const RID &p_joint, PhysicsServer3D::SliderJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_joint_make_cone_twist(const RID &p_joint, const RID &p_body_a,
        const Transform3D &p_local_ref_a, const RID &p_body_b, const Transform3D &p_local_ref_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_CONE_TWIST;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_local_ref_a;
    joint.ref_b = p_local_ref_b;
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_cone_twist_joint_set_param(const RID &p_joint, PhysicsServer3D::ConeTwistJointParam p_param,
        double p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    if (p_param == PhysicsServer3D::CONE_TWIST_JOINT_SWING_SPAN ||
            p_param == PhysicsServer3D::CONE_TWIST_JOINT_TWIST_SPAN) {
        _rebuild_for_joint(p_joint);
    }
}

double AVBDPhysicsServer3D::_cone_twist_joint_get_param(const RID &p_joint,
        PhysicsServer3D::ConeTwistJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_generic_6dof_joint_set_param(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisParam p_param, double p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_param) >= JointData::kMaxParams) {
        return;
    }
    // Rows 0..2 are linear axes, 3..5 angular. The param enum numbers linear and angular
    // families separately; the row offset puts them side by side.
    const bool angular = static_cast<int>(p_param) >= PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT;
    const int param_index = angular
            ? static_cast<int>(p_param) - PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT
            : static_cast<int>(p_param);
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    found->second.params[row][param_index] = p_value;
    if (BodyData *body_b = find_body(found->second.body_b)) {
        rebuild_space(body_b->space);
    } else if (BodyData *body_a = find_body(found->second.body_a)) {
        rebuild_space(body_a->space);
    }
}

double AVBDPhysicsServer3D::_generic_6dof_joint_get_param(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_param) >= JointData::kMaxParams) {
        return 0.0;
    }
    const bool angular = static_cast<int>(p_param) >= PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT;
    const int param_index = angular
            ? static_cast<int>(p_param) - PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT
            : static_cast<int>(p_param);
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    return found->second.params[row][param_index];
}

void AVBDPhysicsServer3D::_generic_6dof_joint_set_flag(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisFlag p_flag, bool p_enable) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_flag) >= JointData::kMaxFlags) {
        return;
    }
    const bool angular = static_cast<int>(p_flag) == PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_ANGULAR_LIMIT ||
            static_cast<int>(p_flag) == g6dof_spring::FLAG_ENABLE_ANGULAR;
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    found->second.flags[row][static_cast<int>(p_flag)] = p_enable;
    if (BodyData *body_b = find_body(found->second.body_b)) {
        rebuild_space(body_b->space);
    } else if (BodyData *body_a = find_body(found->second.body_a)) {
        rebuild_space(body_a->space);
    }
}

bool AVBDPhysicsServer3D::_generic_6dof_joint_get_flag(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisFlag p_flag) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_flag) >= JointData::kMaxFlags) {
        return false;
    }
    const bool angular = static_cast<int>(p_flag) == PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_ANGULAR_LIMIT ||
            static_cast<int>(p_flag) == g6dof_spring::FLAG_ENABLE_ANGULAR;
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    return found->second.flags[row][static_cast<int>(p_flag)];
}

// -----------------------------------------------------------------------------
// Bodies
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_body_create() {
    const uint64_t id = next_id++;
    bodies.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple());
    return make_rid(id);
}

void AVBDPhysicsServer3D::_body_set_space(const RID &p_body, const RID &p_space) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    if (body->space == p_space) {
        return; // pointless, as the built-in server says
    }
    (void)0;
    const RID previous = body->space;
    body->space = p_space;
    // The engine calls this the moment a node enters the tree - before its shapes and
    // physics material arrive. Rebuilding now would freeze an empty body into the
    // solver; deferring to the next step lets the whole setup land first.
    space_rebuild_pending[id_of(p_space)] = true;
    if (previous.is_valid()) {
        space_rebuild_pending[id_of(previous)] = true;
    }
}

RID AVBDPhysicsServer3D::_body_get_space(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? RID() : found->second.space;
}

void AVBDPhysicsServer3D::_body_set_mode(const RID &p_body, PhysicsServer3D::BodyMode p_mode) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->mode = p_mode;
    rebuild_space(body->space);
}

PhysicsServer3D::BodyMode AVBDPhysicsServer3D::_body_get_mode(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? PhysicsServer3D::BODY_MODE_RIGID : found->second.mode;
}

void AVBDPhysicsServer3D::_body_set_state(const RID &p_body, PhysicsServer3D::BodyState p_state, const Variant &p_value) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }

    switch (p_state) {
        case PhysicsServer3D::BODY_STATE_TRANSFORM:
            body->transform = p_value;
            if (body->rigid != nullptr) {
                // A transform set from the scene is the body's new pose right now, so the solver
                // is told to match rather than being left to interpolate towards it.
                const Transform3D &t = body->transform;
                body->rigid->positionLin = avbd::float3{static_cast<float>(t.origin.x),
                        static_cast<float>(-t.origin.z), static_cast<float>(t.origin.y)};
                body->rigid->positionAng = to_sim(t.basis.get_rotation_quaternion());
            }
            break;

        case PhysicsServer3D::BODY_STATE_LINEAR_VELOCITY:
            body->linear_velocity = p_value;
            if (body->rigid != nullptr) {
                const Vector3 &v = body->linear_velocity;
                body->rigid->velocityLin = avbd::float3{static_cast<float>(v.x), static_cast<float>(-v.z),
                        static_cast<float>(v.y)};
            }
            break;

        case PhysicsServer3D::BODY_STATE_ANGULAR_VELOCITY:
            body->angular_velocity = p_value;
            if (body->rigid != nullptr) {
                const Vector3 &w = body->angular_velocity;
                body->rigid->velocityAng = avbd::float3{static_cast<float>(w.x), static_cast<float>(-w.z),
                        static_cast<float>(w.y)};
            }
            break;

        case PhysicsServer3D::BODY_STATE_SLEEPING: {
            body->sleeping = p_value;
            if (body->rigid != nullptr) {
                body->rigid->sleeping = p_value;
                if (p_value) {
                    // Godot zeroes velocity when a body is put to sleep explicitly.
                    body->rigid->velocityLin = avbd::float3{0, 0, 0};
                    body->rigid->velocityAng = avbd::float3{0, 0, 0};
                    body->linear_velocity = Vector3();
                    body->angular_velocity = Vector3();
                }
            }
            break;
        }

        case PhysicsServer3D::BODY_STATE_CAN_SLEEP:
            // Rebuild applies sleep_mode = can_sleep ? SLEEP : NEVER; the setter stores the
            // flag and lets the next rebuild stamp it. Sleeping bodies wake to re-evaluate.
            body->can_sleep = p_value;
            rebuild_space(body->space);
            break;

        default:
            break;
    }
}

Variant AVBDPhysicsServer3D::_body_get_state(const RID &p_body, PhysicsServer3D::BodyState p_state) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return Variant();
    }
    const BodyData &body = found->second;

    switch (p_state) {
        case PhysicsServer3D::BODY_STATE_TRANSFORM:
            return body.transform;
        case PhysicsServer3D::BODY_STATE_LINEAR_VELOCITY:
            return body.linear_velocity;
        case PhysicsServer3D::BODY_STATE_ANGULAR_VELOCITY:
            return body.angular_velocity;
        case PhysicsServer3D::BODY_STATE_SLEEPING:
            return body.rigid != nullptr ? body.rigid->sleeping : body.sleeping;
        case PhysicsServer3D::BODY_STATE_CAN_SLEEP:
            return body.can_sleep;
        default:
            return false;
    }
}

void AVBDPhysicsServer3D::_body_set_state_sync_callback(const RID &p_body, const Callable &p_callable) {
    if (BodyData *body = find_body(p_body)) {
        body->sync = p_callable;
    }
}

void AVBDPhysicsServer3D::_body_attach_object_instance_id(const RID &p_body, uint64_t p_id) {
    if (BodyData *body = find_body(p_body)) {
        body->instance_id = p_id;
    }
}

uint64_t AVBDPhysicsServer3D::_body_get_object_instance_id(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : found->second.instance_id;
}

void AVBDPhysicsServer3D::_body_set_ray_pickable(const RID &p_body, bool p_enable) {
    // Space queries do not read this yet; storing keeps the getter honest.
    if (BodyData *body = find_body(p_body)) {
        body->ray_pickable = p_enable;
    }
}

// --- DirectBodyState accessors -------------------------------------------------
//
// The state object is the engine's view of one body during (or after) a step; these
// helpers do the BodyData lookups it cannot do itself. All of them tolerate an unknown
// RID: the engine may ask about a body this server never saw.

void AVBDPhysicsServer3D::body_state_add_constant_force(const RID &p_body, const Vector3 &p_force,
        const Vector3 &p_position) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force += p_force;
        if (p_position.length_squared() > 0.0) {
            body->constant_torque += p_position.cross(p_force);
        }
    }
}

void AVBDPhysicsServer3D::body_state_add_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque += p_torque;
    }
}

void AVBDPhysicsServer3D::body_state_set_constant_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force = p_force;
    }
}

Vector3 AVBDPhysicsServer3D::body_state_get_constant_force(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_force;
}

void AVBDPhysicsServer3D::body_state_set_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque = p_torque;
    }
}

Vector3 AVBDPhysicsServer3D::body_state_get_constant_torque(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_torque;
}

int32_t AVBDPhysicsServer3D::body_state_contact_count(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : static_cast<int32_t>(found->second.frame_contacts.size());
}

Vector3 AVBDPhysicsServer3D::body_state_contact_position(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return Vector3();
    }
    return found->second.frame_contacts[p_index].position;
}

RID AVBDPhysicsServer3D::body_state_contact_collider(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return RID();
    }
    return make_rid(found->second.frame_contacts[p_index].collider_body_id);
}

uint64_t AVBDPhysicsServer3D::body_state_contact_collider_id(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return 0;
    }
    const BodyData *other = find_body(make_rid(found->second.frame_contacts[p_index].collider_body_id));
    return other != nullptr ? other->instance_id : 0;
}

Vector3 AVBDPhysicsServer3D::body_state_contact_collider_velocity(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return Vector3();
    }
    const BodyData *other = find_body(make_rid(found->second.frame_contacts[p_index].collider_body_id));
    return other != nullptr ? other->linear_velocity : Vector3();
}

double AVBDPhysicsServer3D::body_state_step() const {
    // Every space steps with the engine's fixed timestep; report the common value.
    for (const auto &[id, space] : spaces) {
        (void)id;
        return space.solver.dt;
    }
    return 0.0;
}

PhysicsDirectSpaceState3D *AVBDPhysicsServer3D::body_state_space_state(const RID &p_body) {
    const BodyData *body = find_body(p_body);
    if (body == nullptr || body->space.is_valid() == false) {
        return nullptr;
    }
    return _space_get_direct_state(body->space);
}

PhysicsDirectBodyState3D *AVBDPhysicsServer3D::_body_get_direct_state(const RID &p_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return nullptr;
    }
    if (body_state == nullptr) {
        body_state = memnew(AVBDDirectBodyState3D);
    }
    body_state->bind(this, body->rigid, p_body);
    return body_state;
}

// --- query helpers -------------------------------------------------------------

std::vector<AVBDPhysicsServer3D::QueryCandidate> AVBDPhysicsServer3D::query_candidates(const RID &p_space,
        uint32_t p_mask) const {
    std::vector<QueryCandidate> out;
    const auto found = spaces.find(id_of(p_space));
    if (found == spaces.end()) {
        return out;
    }
    const uint64_t space_id = id_of(p_space);
    for (const auto &[id, body] : bodies) {
        if (body.rigid == nullptr || id_of(body.space) != space_id) {
            continue;
        }
        if ((body.rigid->collisionLayer & p_mask) == 0) {
            continue;
        }
        out.push_back({body.rigid, id});
    }
    return out;
}

avbd::Rigid *AVBDPhysicsServer3D::solver_pick(const RID &p_space, avbd::float3 p_origin, avbd::float3 p_dir,
        uint32_t p_mask, avbd::float3 &r_local) {
    SpaceData *space = find_space(p_space);
    if (space == nullptr) {
        return nullptr;
    }
    return space->solver.pick(p_origin, p_dir, r_local, p_mask);
}

uint64_t AVBDPhysicsServer3D::body_id_of(const RID &p_space, const avbd::Rigid *p_rigid) const {
    const uint64_t space_id = id_of(p_space);
    for (const auto &[id, body] : bodies) {
        if (body.rigid == p_rigid && id_of(body.space) == space_id) {
            return id;
        }
    }
    return 0;
}

RID AVBDPhysicsServer3D::solver_pick_rid(const RID &p_space, const avbd::Rigid *p_rigid) {
    const uint64_t id = body_id_of(p_space, p_rigid);
    return id != 0 ? make_rid(id) : RID();
}

uint64_t AVBDPhysicsServer3D::body_instance_id(const RID &p_space, const avbd::Rigid *p_rigid) {
    const BodyData *body = find_body(make_rid(body_id_of(p_space, p_rigid)));
    return body != nullptr ? body->instance_id : 0;
}

Vector3 AVBDPhysicsServer3D::body_velocity_of(const RID &p_space, const avbd::Rigid *p_rigid) {
    const BodyData *body = find_body(make_rid(body_id_of(p_space, p_rigid)));
    return body != nullptr ? body->linear_velocity : Vector3();
}

const AVBDPhysicsServer3D::ShapeData *AVBDPhysicsServer3D::shape_data(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    return found == shapes.end() ? nullptr : &found->second;
}

// -----------------------------------------------------------------------------
// Body parameters and forces
// -----------------------------------------------------------------------------

void AVBDPhysicsServer3D::_body_set_param(const RID &p_body, PhysicsServer3D::BodyParameter p_param,
        const Variant &p_value) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }

    switch (p_param) {
        case PhysicsServer3D::BODY_PARAM_MASS:
            body->mass = p_value;
            rebuild_space(body->space);
            break;
        case PhysicsServer3D::BODY_PARAM_FRICTION:
            body->friction = p_value;
            if (body->rigid != nullptr) {
                body->rigid->friction = body->friction;
            }
            break;
        case PhysicsServer3D::BODY_PARAM_BOUNCE:
            // The AVBD contact model is friction-only: a bounce value is remembered for
            // _body_get_param but no restitution is simulated.
            body->bounce = p_value;
            break;
        case PhysicsServer3D::BODY_PARAM_GRAVITY_SCALE: {
            body->gravity_scale = p_value;
            if (body->rigid != nullptr && body->space.is_valid()) {
                // Re-derive the per-body gravity from the space's magnitude.
                if (SpaceData *space = find_space(body->space)) {
                    body->rigid->gravity = space->solver.gravity * static_cast<float>(body->gravity_scale);
                }
            }
            break;
        }
        case PhysicsServer3D::BODY_PARAM_LINEAR_DAMP:
            body->linear_damp = p_value;
            break; // damping is not simulated (documented)
        case PhysicsServer3D::BODY_PARAM_ANGULAR_DAMP:
            body->angular_damp = p_value;
            break; // damping is not simulated (documented)
        case PhysicsServer3D::BODY_PARAM_LINEAR_DAMP_MODE:
        case PhysicsServer3D::BODY_PARAM_ANGULAR_DAMP_MODE:
            break; // damp modes are meaningless without damping
        case PhysicsServer3D::BODY_PARAM_INERTIA:
        case PhysicsServer3D::BODY_PARAM_CENTER_OF_MASS:
            break; // inertia follows the shape; the centre of mass is the body origin
        default:
            break;
    }
}

Variant AVBDPhysicsServer3D::_body_get_param(const RID &p_body, PhysicsServer3D::BodyParameter p_param) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return Variant();
    }
    const BodyData &body = found->second;

    switch (p_param) {
        case PhysicsServer3D::BODY_PARAM_MASS: {
            if (body.rigid != nullptr && body.rigid->mass > 0.0f) {
                return body.rigid->mass; // the solver's value, from density * volume
            }
            return body.mass;
        }
        case PhysicsServer3D::BODY_PARAM_FRICTION:
            return body.rigid != nullptr ? body.rigid->friction : body.friction;
        case PhysicsServer3D::BODY_PARAM_BOUNCE:
            return body.bounce;
        case PhysicsServer3D::BODY_PARAM_GRAVITY_SCALE:
            return body.gravity_scale;
        case PhysicsServer3D::BODY_PARAM_LINEAR_DAMP:
            return body.linear_damp;
        case PhysicsServer3D::BODY_PARAM_ANGULAR_DAMP:
            return body.angular_damp;
        default:
            return 0.0;
    }
}

// Impulses and forces arrive in Godot's world frame; the solver lives in its own, so everything
// below goes through to_sim. The math mirrors the node layer's apply_impulse/apply_torque_impulse.

void AVBDPhysicsServer3D::_body_apply_central_impulse(const RID &p_body, const Vector3 &p_impulse) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || body->rigid == nullptr || body->rigid->mass <= 0.0f) {
        return;
    }
    body->rigid->velocityLin += to_sim(p_impulse) / body->rigid->mass;
}

void AVBDPhysicsServer3D::_body_apply_impulse(const RID &p_body, const Vector3 &p_impulse,
        const Vector3 &p_position) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || body->rigid == nullptr || body->rigid->mass <= 0.0f) {
        return;
    }
    avbd::Rigid &rigid = *body->rigid;
    const avbd::float3 impulse = to_sim(p_impulse);
    rigid.velocityLin += impulse / rigid.mass;

    if (p_position.length_squared() > 0.0) {
        const avbd::float3 torque = avbd::cross(to_sim(p_position), impulse);
        const avbd::float3 local = avbd::rotate(avbd::conjugate(rigid.positionAng), torque);
        const avbd::float3 delta{local.x / rigid.moment.x, local.y / rigid.moment.y,
                local.z / rigid.moment.z};
        rigid.velocityAng += avbd::rotate(rigid.positionAng, delta);
    }
}

void AVBDPhysicsServer3D::_body_apply_torque_impulse(const RID &p_body, const Vector3 &p_impulse) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || body->rigid == nullptr || body->rigid->mass <= 0.0f) {
        return;
    }
    avbd::Rigid &rigid = *body->rigid;
    const avbd::float3 local = avbd::rotate(avbd::conjugate(rigid.positionAng), to_sim(p_impulse));
    const avbd::float3 delta{local.x / rigid.moment.x, local.y / rigid.moment.y,
            local.z / rigid.moment.z};
    rigid.velocityAng += avbd::rotate(rigid.positionAng, delta);
}

void AVBDPhysicsServer3D::_body_apply_central_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->frame_force += p_force;
    }
}

void AVBDPhysicsServer3D::_body_apply_force(const RID &p_body, const Vector3 &p_force,
        const Vector3 &p_position) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->frame_force += p_force;
    if (p_position.length_squared() > 0.0) {
        body->frame_torque += p_position.cross(p_force);
    }
}

void AVBDPhysicsServer3D::_body_apply_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->frame_torque += p_torque;
    }
}

void AVBDPhysicsServer3D::_body_add_constant_central_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force = p_force;
    }
}

void AVBDPhysicsServer3D::_body_add_constant_force(const RID &p_body, const Vector3 &p_force,
        const Vector3 &p_position) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->constant_force += p_force;
    if (p_position.length_squared() > 0.0) {
        body->constant_torque += p_position.cross(p_force);
    }
}

void AVBDPhysicsServer3D::_body_add_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque += p_torque;
    }
}

// Godot BodyAxis bits (linear X/Y/Z = 1/2/4, angular X/Y/Z = 8/16/32) to the solver's
// sim-axis bits (0..2 = sim x/y/z). Godot Y is sim Z, Godot Z is sim -Y; the lock is a
// directionless bit so the sign does not matter, only the permutation.
void AVBDPhysicsServer3D::_body_set_axis_lock(const RID &p_body, PhysicsServer3D::BodyAxis p_axis, bool p_lock) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    const uint32_t bit = static_cast<uint32_t>(p_axis);
    if (p_lock) {
        body->axis_locks |= bit;
    } else {
        body->axis_locks &= ~bit;
    }
    if (body->rigid != nullptr) {
        uint8_t linear = 0;
        uint8_t angular = 0;
        axis_locks_to_sim(body->axis_locks, linear, angular);
        body->rigid->axisLockLinear = linear;
        body->rigid->axisLockAngular = angular;
    }
}

bool AVBDPhysicsServer3D::_body_is_axis_locked(const RID &p_body, PhysicsServer3D::BodyAxis p_axis) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return false;
    }
    return (found->second.axis_locks & static_cast<uint32_t>(p_axis)) != 0;
}

void AVBDPhysicsServer3D::_body_add_collision_exception(const RID &p_body, const RID &p_excepted_body) {
    BodyData *body = find_body(p_body);
    const uint64_t other = id_of(p_excepted_body);
    if (body == nullptr || other == id_of(p_body)) {
        return;
    }
    if (std::find(body->exceptions.begin(), body->exceptions.end(), other) == body->exceptions.end()) {
        body->exceptions.push_back(other);
    }
    if (body->space.is_valid()) {
        rebuild_space(body->space);
    }
}

void AVBDPhysicsServer3D::_body_remove_collision_exception(const RID &p_body, const RID &p_excepted_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    std::vector<uint64_t> &list = body->exceptions;
    list.erase(std::remove(list.begin(), list.end(), id_of(p_excepted_body)), list.end());
    if (body->space.is_valid()) {
        rebuild_space(body->space);
    }
}

TypedArray<RID> AVBDPhysicsServer3D::_body_get_collision_exceptions(const RID &p_body) const {
    TypedArray<RID> out;
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return out;
    }
    for (const uint64_t other_id : found->second.exceptions) {
        out.push_back(make_rid(other_id));
    }
    return out;
}

void AVBDPhysicsServer3D::_body_set_max_contacts_reported(const RID &p_body, int32_t p_amount) {
    if (BodyData *body = find_body(p_body)) {
        body->max_contacts_reported = p_amount;
    }
}

int32_t AVBDPhysicsServer3D::_body_get_max_contacts_reported(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : found->second.max_contacts_reported;
}

void AVBDPhysicsServer3D::_body_set_contacts_reported_depth_threshold(const RID &p_body, double p_threshold) {
    if (BodyData *body = find_body(p_body)) {
        body->contacts_depth_threshold = p_threshold; // recorded; not applied (contact depths come from the solver)
    }
}

double AVBDPhysicsServer3D::_body_get_contacts_reported_depth_threshold(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0.0 : found->second.contacts_depth_threshold;
}

void AVBDPhysicsServer3D::_body_set_omit_force_integration(const RID &p_body, bool p_enable) {
    if (BodyData *body = find_body(p_body)) {
        body->omit_force_integration = p_enable;
    }
}

bool AVBDPhysicsServer3D::_body_is_omitting_force_integration(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? false : found->second.omit_force_integration;
}

void AVBDPhysicsServer3D::_body_set_force_integration_callback(const RID &p_body, const Callable &p_callable,
        const Variant &p_userdata) {
    if (BodyData *body = find_body(p_body)) {
        body->force_integration_callback = p_callable;
        body->force_integration_userdata = p_userdata;
    }
}

void AVBDPhysicsServer3D::_body_set_constant_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force = p_force;
    }
}

Vector3 AVBDPhysicsServer3D::_body_get_constant_force(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_force;
}

void AVBDPhysicsServer3D::_body_set_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque = p_torque;
    }
}

Vector3 AVBDPhysicsServer3D::_body_get_constant_torque(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_torque;
}

void AVBDPhysicsServer3D::_body_set_axis_velocity(const RID &p_body, const Vector3 &p_axis_velocity) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    // Godot semantics: each nonzero component replaces the velocity along that axis.
    Vector3 velocity = body->linear_velocity;
    if (p_axis_velocity.x != 0.0) {
        velocity.x = p_axis_velocity.x;
    }
    if (p_axis_velocity.y != 0.0) {
        velocity.y = p_axis_velocity.y;
    }
    if (p_axis_velocity.z != 0.0) {
        velocity.z = p_axis_velocity.z;
    }
    body->linear_velocity = velocity;
    if (body->rigid != nullptr) {
        body->rigid->velocityLin = to_sim(velocity);
    }
}

void AVBDPhysicsServer3D::_body_set_enable_continuous_collision_detection(const RID &p_body, bool p_enable) {
    if (BodyData *body = find_body(p_body)) {
        body->ccd = p_enable; // recorded; CCD is not simulated (documented)
    }
}

bool AVBDPhysicsServer3D::_body_is_continuous_collision_detection_enabled(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? false : found->second.ccd;
}

void AVBDPhysicsServer3D::_body_set_collision_layer(const RID &p_body, uint32_t p_layer) {
    if (BodyData *body = find_body(p_body)) {
        body->collision_layer = p_layer;
        if (body->rigid != nullptr) {
            body->rigid->collisionLayer = p_layer;
        }
    }
}

uint32_t AVBDPhysicsServer3D::_body_get_collision_layer(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 1 : found->second.collision_layer;
}

void AVBDPhysicsServer3D::_body_set_collision_mask(const RID &p_body, uint32_t p_mask) {
    if (BodyData *body = find_body(p_body)) {
        body->collision_mask = p_mask;
        if (body->rigid != nullptr) {
            body->rigid->collisionMask = p_mask;
        }
    }
}

uint32_t AVBDPhysicsServer3D::_body_get_collision_mask(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 1 : found->second.collision_mask;
}

void AVBDPhysicsServer3D::_body_set_collision_priority(const RID &p_body, double p_priority) {
    if (BodyData *body = find_body(p_body)) {
        body->collision_priority = p_priority; // recorded; solver has no priority weighting
    }
}

double AVBDPhysicsServer3D::_body_get_collision_priority(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 1.0 : found->second.collision_priority;
}

void AVBDPhysicsServer3D::_body_set_user_flags(const RID &p_body, uint32_t p_flags) {
    if (BodyData *body = find_body(p_body)) {
        body->user_flags = p_flags;
    }
}

uint32_t AVBDPhysicsServer3D::_body_get_user_flags(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : found->second.user_flags;
}

void AVBDPhysicsServer3D::_body_reset_mass_properties(const RID &p_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    // Godot's contract for this call is "restore the default inertia and centre of mass to
    // cancel custom values". Both are always derived from the shape here, so the rebuild is
    // the whole of the reset - it must not touch the engine's mass.
    rebuild_space(body->space);
}

void AVBDPhysicsServer3D::_body_add_shape(const RID &p_body, const RID &p_shape, const Transform3D &p_transform,
        bool) {
    BodyData *body = find_body(p_body);
    ShapeData *shape = find_shape(p_shape);
    if (body == nullptr || shape == nullptr) {
        return;
    }
    // Keep the RID: the engine may set the shape resource's data after attaching it, and
    // the deferred rebuild refreshes the snapshot from this RID.
    body->shape_rids.push_back(p_shape);
    body->shapes.push_back(*shape);
    body->shape_transforms.push_back(p_transform);
    space_rebuild_pending[id_of(body->space)] = true;
}

void AVBDPhysicsServer3D::_body_set_shape(const RID &p_body, int32_t p_shape_idx, const RID &p_shape) {
    BodyData *body = find_body(p_body);
    ShapeData *shape = find_shape(p_shape);
    if (body == nullptr || shape == nullptr || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= body->shapes.size()) {
        return;
    }
    body->shape_rids[p_shape_idx] = p_shape;
    body->shapes[p_shape_idx] = *shape;
    space_rebuild_pending[id_of(body->space)] = true;
}

void AVBDPhysicsServer3D::_body_set_shape_transform(const RID &p_body, int32_t p_shape_idx,
        const Transform3D &p_transform) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= body->shape_transforms.size()) {
        return;
    }
    body->shape_transforms[p_shape_idx] = p_transform;
}

void AVBDPhysicsServer3D::_body_set_shape_disabled(const RID &p_body, int32_t p_shape_idx, bool p_disabled) {
    // A disabled shape contributes nothing; AVBD has one shape per body, so a disabled one leaves
    // the body without collision rather than with a smaller set.
    BodyData *body = find_body(p_body);
    if (body == nullptr || p_shape_idx < 0) {
        return;
    }
    if (p_disabled) {
        body->shapes.clear();
        body->shape_transforms.clear();
        body->shape_rids.clear();
        space_rebuild_pending[id_of(body->space)] = true;
    }
}

int32_t AVBDPhysicsServer3D::_body_get_shape_count(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : static_cast<int32_t>(found->second.shapes.size());
}

RID AVBDPhysicsServer3D::_body_get_shape(const RID &p_body, int32_t p_shape_idx) const {
    // The engine keeps the shape RIDs; this server only records the parameters it was given, so
    // there is no handle to give back. Returning an invalid RID is honest - the caller is asking
    // for something this implementation does not track.
    (void)p_body;
    (void)p_shape_idx;
    return RID();
}

Transform3D AVBDPhysicsServer3D::_body_get_shape_transform(const RID &p_body, int32_t p_shape_idx) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= found->second.shape_transforms.size()) {
        return Transform3D();
    }
    return found->second.shape_transforms[p_shape_idx];
}

void AVBDPhysicsServer3D::_body_remove_shape(const RID &p_body, int32_t p_shape_idx) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= body->shapes.size()) {
        return;
    }
    body->shapes.erase(body->shapes.begin() + p_shape_idx);
    body->shape_transforms.erase(body->shape_transforms.begin() + p_shape_idx);
    body->shape_rids.erase(body->shape_rids.begin() + p_shape_idx);
    space_rebuild_pending[id_of(body->space)] = true;
}

void AVBDPhysicsServer3D::_body_clear_shapes(const RID &p_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->shapes.clear();
    body->shape_transforms.clear();
    body->shape_rids.clear();
    space_rebuild_pending[id_of(body->space)] = true;
}

// -----------------------------------------------------------------------------
// Areas
//
// Every storage virtual keeps the engine's description in an AreaData record, and
// _monitor_areas diffs each area against the bodies of its space every step, calling
// the monitor callbacks with Godot's four-argument contract (verified against 4.7
// scene/3d/physics/area_3d.cpp): (status, body_rid, instance_id, body_shape, area_shape)
// - status 0 = entered (AREA_BODY_ADDED), 1 = exited (AREA_BODY_REMOVED).
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_area_create() {
    const uint64_t id = next_id++;
    areas.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple());
    return make_rid(id);
}

AVBDPhysicsServer3D::AreaData *AVBDPhysicsServer3D::find_area(const RID &p_rid) {
    const auto found = areas.find(id_of(p_rid));
    return found == areas.end() ? nullptr : &found->second;
}

void AVBDPhysicsServer3D::_area_set_space(const RID &p_area, const RID &p_space) {
    if (AreaData *area = find_area(p_area)) {
        area->space = p_space;
        area->prev_inside_bodies.clear();
        area->prev_inside_areas.clear();
    }
}

RID AVBDPhysicsServer3D::_area_get_space(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? RID() : found->second.space;
}

void AVBDPhysicsServer3D::_area_set_transform(const RID &p_area, const Transform3D &p_transform) {
    if (AreaData *area = find_area(p_area)) {
        area->transform = p_transform;
    }
}

Transform3D AVBDPhysicsServer3D::_area_get_transform(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? Transform3D() : found->second.transform;
}

void AVBDPhysicsServer3D::_area_add_shape(const RID &p_area, const RID &p_shape, const Transform3D &p_transform,
        bool p_disabled) {
    AreaData *area = find_area(p_area);
    ShapeData *shape = find_shape(p_shape);
    if (area == nullptr || shape == nullptr) {
        return;
    }
    area->shapes.push_back(*shape);
    area->shape_transforms.push_back(p_transform);
    area->shape_disabled.push_back(p_disabled);
}

void AVBDPhysicsServer3D::_area_set_shape(const RID &p_area, int32_t p_shape_idx, const RID &p_shape) {
    AreaData *area = find_area(p_area);
    ShapeData *shape = find_shape(p_shape);
    if (area == nullptr || shape == nullptr || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= area->shapes.size()) {
        return;
    }
    area->shapes[p_shape_idx] = *shape;
}

void AVBDPhysicsServer3D::_area_set_shape_transform(const RID &p_area, int32_t p_shape_idx,
        const Transform3D &p_transform) {
    AreaData *area = find_area(p_area);
    if (area == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= area->shape_transforms.size()) {
        return;
    }
    area->shape_transforms[p_shape_idx] = p_transform;
}

void AVBDPhysicsServer3D::_area_set_shape_disabled(const RID &p_area, int32_t p_shape_idx, bool p_disabled) {
    AreaData *area = find_area(p_area);
    if (area == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= area->shape_disabled.size()) {
        return;
    }
    area->shape_disabled[p_shape_idx] = p_disabled;
}

int32_t AVBDPhysicsServer3D::_area_get_shape_count(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 0 : static_cast<int32_t>(found->second.shapes.size());
}

RID AVBDPhysicsServer3D::_area_get_shape(const RID &, int32_t) const {
    // As with bodies: the engine owns the shape RIDs; this server only keeps parameters.
    return RID();
}

Transform3D AVBDPhysicsServer3D::_area_get_shape_transform(const RID &p_area, int32_t p_shape_idx) const {
    const auto found = areas.find(id_of(p_area));
    if (found == areas.end() || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= found->second.shape_transforms.size()) {
        return Transform3D();
    }
    return found->second.shape_transforms[p_shape_idx];
}

void AVBDPhysicsServer3D::_area_remove_shape(const RID &p_area, int32_t p_shape_idx) {
    AreaData *area = find_area(p_area);
    if (area == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= area->shapes.size()) {
        return;
    }
    area->shapes.erase(area->shapes.begin() + p_shape_idx);
    area->shape_transforms.erase(area->shape_transforms.begin() + p_shape_idx);
    area->shape_disabled.erase(area->shape_disabled.begin() + p_shape_idx);
}

void AVBDPhysicsServer3D::_area_clear_shapes(const RID &p_area) {
    if (AreaData *area = find_area(p_area)) {
        area->shapes.clear();
        area->shape_transforms.clear();
        area->shape_disabled.clear();
    }
}

void AVBDPhysicsServer3D::_area_set_param(const RID &p_area, PhysicsServer3D::AreaParameter p_param,
        const Variant &p_value) {
    if (AreaData *area = find_area(p_area)) {
        // Recorded so the getter replays it; gravity overrides and the other area physics
        // parameters do not act on the solver this round (documented).
        area->params[p_param] = p_value;
    }
}

Variant AVBDPhysicsServer3D::_area_get_param(const RID &p_area, PhysicsServer3D::AreaParameter p_param) const {
    const auto found = areas.find(id_of(p_area));
    if (found == areas.end()) {
        return Variant();
    }
    const auto param = found->second.params.find(p_param);
    return param == found->second.params.end() ? Variant() : param->second;
}

void AVBDPhysicsServer3D::_area_set_collision_layer(const RID &p_area, uint32_t p_layer) {
    if (AreaData *area = find_area(p_area)) {
        area->collision_layer = p_layer;
    }
}

uint32_t AVBDPhysicsServer3D::_area_get_collision_layer(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 1 : found->second.collision_layer;
}

void AVBDPhysicsServer3D::_area_set_collision_mask(const RID &p_area, uint32_t p_mask) {
    if (AreaData *area = find_area(p_area)) {
        area->collision_mask = p_mask;
    }
}

uint32_t AVBDPhysicsServer3D::_area_get_collision_mask(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 1 : found->second.collision_mask;
}

void AVBDPhysicsServer3D::_area_set_monitorable(const RID &p_area, bool p_monitorable) {
    if (AreaData *area = find_area(p_area)) {
        area->monitorable = p_monitorable;
    }
}

void AVBDPhysicsServer3D::_area_set_ray_pickable(const RID &p_area, bool p_enable) {
    if (AreaData *area = find_area(p_area)) {
        area->ray_pickable = p_enable; // space queries do not hit areas this round
    }
}

void AVBDPhysicsServer3D::_area_attach_object_instance_id(const RID &p_area, uint64_t p_id) {
    if (AreaData *area = find_area(p_area)) {
        area->instance_id = p_id;
    }
}

uint64_t AVBDPhysicsServer3D::_area_get_object_instance_id(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 0 : found->second.instance_id;
}

void AVBDPhysicsServer3D::_area_set_monitor_callback(const RID &p_area, const Callable &p_callback) {
    if (AreaData *area = find_area(p_area)) {
        area->monitor_callback = p_callback;
        area->monitor_callback_set = p_callback.is_valid();
    }
}

void AVBDPhysicsServer3D::_area_set_area_monitor_callback(const RID &p_area, const Callable &p_callback) {
    if (AreaData *area = find_area(p_area)) {
        area->area_monitor_callback = p_callback;
        area->area_monitor_callback_set = p_callback.is_valid();
    }
}

// One step of area monitoring: diff every area's inside-set against the previous one and
// fire the callbacks on the transitions. Godot's body-shape/area-shape arguments are the
// first overlapping pair (this server keeps one solver shape per body, so 0 on the body
// side and the area's first enabled shape on its side).
void AVBDPhysicsServer3D::_monitor_areas() {
    for (auto &[area_id, area] : areas) {
        (void)area_id;
        if (!area.space.is_valid()) {
            continue;
        }
        const SpaceData *space = find_space(area.space);
        if (space == nullptr) {
            continue;
        }

        // The area's shapes in sim space, once per area.
        std::vector<avbd::Shape> area_shapes;
        for (size_t i = 0; i < area.shapes.size(); i++) {
            if (area.shape_disabled[i]) {
                continue;
            }
            const Transform3D posed = area.transform * area.shape_transforms[i];
            area_shapes.push_back([this, &area, i, posed]() {
                avbd::Shape s;
                const ShapeData &data = area.shapes[i];
                switch (data.type) {
                    case PhysicsServer3D::SHAPE_SPHERE:
                        s.type = avbd::ShapeType::Sphere;
                        s.radius = static_cast<float>(data.radius);
                        break;
                    case PhysicsServer3D::SHAPE_CYLINDER:
                        s.type = avbd::ShapeType::Cylinder;
                        s.radius = static_cast<float>(data.radius);
                        s.halfHeight = static_cast<float>(data.height) * 0.5f;
                        break;
                    default:
                        s.type = avbd::ShapeType::Box;
                        s.half = to_sim_extents(data.extents) * 0.5f;
                        break;
                }
                s.center = to_sim(posed.origin);
                s.rotation = to_sim(posed.basis.get_rotation_quaternion());
                s.axis = avbd::rotate(s.rotation, avbd::float3{0, 0, 1});
                return s;
            }());
        }
        if (area_shapes.empty()) {
            continue;
        }

        // --- bodies: the area reports what enters its mask ---
        std::vector<uint64_t> inside;
        if (area.monitor_callback_set || area.monitor_callback.is_valid()) {
            for (const auto &[body_id, body] : bodies) {
                if (body.rigid == nullptr || id_of(body.space) != id_of(area.space)) {
                    continue;
                }
                // Godot's area test: the body's layer must meet the area's mask.
                if ((body.rigid->collisionLayer & area.collision_mask) == 0) {
                    continue;
                }
                const avbd::Shape body_shape = [r = body.rigid]() {
                    avbd::Shape s;
                    s.type = r->shape;
                    s.center = r->positionLin;
                    s.rotation = r->positionAng;
                    s.half = r->size * 0.5f;
                    s.radius = r->size.x;
                    s.halfHeight = r->size.z * 0.5f;
                    s.axis = avbd::rotate(r->positionAng, avbd::float3{0, 0, 1});
                    return s;
                }();
                for (const avbd::Shape &shape : area_shapes) {
                    avbd::Manifold::Contact contact{};
                    avbd::float3x3 basis{};
                    if (avbd::collideShapes(shape, body_shape, &contact, basis) > 0) {
                        inside.push_back(body_id);
                        break;
                    }
                }
            }

            // Diff against the previous step and fire.
            auto was_inside = [&area](uint64_t id) -> bool {
                return std::find(area.prev_inside_bodies.begin(), area.prev_inside_bodies.end(), id)
                        != area.prev_inside_bodies.end();
            };
            for (const uint64_t body_id : inside) {
                if (was_inside(body_id) || !area.monitor_callback.is_valid()) {
                    continue;
                }
                const BodyData *body = find_body(make_rid(body_id));
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_ADDED));
                args.push_back(make_rid(body_id));
                args.push_back(body != nullptr ? body->instance_id : 0);
                args.push_back(0); // body shape index (one shape per solver body)
                args.push_back(0); // area shape index (first shape)
                area.monitor_callback.callv(args);
            }
            for (const uint64_t body_id : area.prev_inside_bodies) {
                if (std::find(inside.begin(), inside.end(), body_id) != inside.end()) {
                    continue;
                }
                if (!area.monitor_callback.is_valid()) {
                    break;
                }
                const BodyData *body = find_body(make_rid(body_id));
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_REMOVED));
                args.push_back(make_rid(body_id));
                args.push_back(body != nullptr ? body->instance_id : 0);
                args.push_back(0);
                args.push_back(0);
                area.monitor_callback.callv(args);
            }
            area.prev_inside_bodies = inside;
        } else {
            area.prev_inside_bodies.clear();
        }

        // --- areas monitoring areas (only when monitorable) ---
        if (area.area_monitor_callback.is_valid()) {
            std::vector<uint64_t> areas_inside;
            if (area.monitorable || true) { // area×area detection is symmetric and cheap
                for (auto &[other_id, other] : areas) {
                    if (other_id == area_id || !other.space.is_valid() || id_of(other.space) != id_of(area.space)) {
                        continue;
                    }
                    if (!other.monitorable) {
                        continue; // Godot: only monitorable areas are seen by other areas
                    }
                    // The other area's layer must meet this area's mask.
                    if ((other.collision_layer & area.collision_mask) == 0) {
                        continue;
                    }
                    bool overlaps = false;
                    for (size_t i = 0; i < other.shapes.size() && !overlaps; i++) {
                        if (other.shape_disabled[i]) {
                            continue;
                        }
                        const Transform3D posed = other.transform * other.shape_transforms[i];
                        const ShapeData &data = other.shapes[i];
                        avbd::Shape other_shape;
                        switch (data.type) {
                            case PhysicsServer3D::SHAPE_SPHERE:
                                other_shape.type = avbd::ShapeType::Sphere;
                                other_shape.radius = static_cast<float>(data.radius);
                                break;
                            case PhysicsServer3D::SHAPE_CYLINDER:
                                other_shape.type = avbd::ShapeType::Cylinder;
                                other_shape.radius = static_cast<float>(data.radius);
                                other_shape.halfHeight = static_cast<float>(data.height) * 0.5f;
                                break;
                            default:
                                other_shape.type = avbd::ShapeType::Box;
                                other_shape.half = to_sim_extents(data.extents) * 0.5f;
                                break;
                        }
                        other_shape.center = to_sim(posed.origin);
                        other_shape.rotation = to_sim(posed.basis.get_rotation_quaternion());
                        other_shape.axis = avbd::rotate(other_shape.rotation, avbd::float3{0, 0, 1});

                        for (const avbd::Shape &shape : area_shapes) {
                            avbd::Manifold::Contact contact{};
                            avbd::float3x3 basis{};
                            if (avbd::collideShapes(shape, other_shape, &contact, basis) > 0) {
                                overlaps = true;
                                break;
                            }
                        }
                    }
                    if (overlaps) {
                        areas_inside.push_back(other_id);
                    }
                }
            }

            auto area_was_inside = [&areas_inside, &area](uint64_t id) -> bool {
                return std::find(area.prev_inside_areas.begin(), area.prev_inside_areas.end(), id)
                        != area.prev_inside_areas.end();
            };
            for (const uint64_t other_id : areas_inside) {
                if (area_was_inside(other_id)) {
                    continue;
                }
                const auto other_found = areas.find(other_id);
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_ADDED));
                args.push_back(make_rid(other_id));
                args.push_back(other_found != areas.end() ? other_found->second.instance_id : 0);
                args.push_back(0);
                args.push_back(0);
                area.area_monitor_callback.callv(args);
            }
            for (const uint64_t other_id : area.prev_inside_areas) {
                if (std::find(areas_inside.begin(), areas_inside.end(), other_id) != areas_inside.end()) {
                    continue;
                }
                const auto other_found = areas.find(other_id);
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_REMOVED));
                args.push_back(make_rid(other_id));
                args.push_back(other_found != areas.end() ? other_found->second.instance_id : 0);
                args.push_back(0);
                args.push_back(0);
                area.area_monitor_callback.callv(args);
            }
            area.prev_inside_areas = areas_inside;
        } else {
            area.prev_inside_areas.clear();
        }
    }
}

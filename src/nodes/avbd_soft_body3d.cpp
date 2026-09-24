/*
 * AVBDSoftBody3D - a lattice of jointed boxes that behaves like a soft body.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#include "nodes/avbd_soft_body3d.hpp"

#include <godot_cpp/classes/box_mesh.hpp>
#include <godot_cpp/classes/multi_mesh.hpp>
#include <algorithm>
#include <cmath>
#include <ranges>

#include <godot_cpp/core/class_db.hpp>

#include "nodes/godot_convert.hpp"

using namespace godot;
using namespace avbd_godot;

namespace {

// Godot axis -> solver axis, for the +90 degree rotation about X used everywhere
// in this module (Godot x -> sim x, Godot y -> sim z, Godot z -> sim y).
constexpr int SIM_AXIS[3] = {0, 2, 1};

// Offset of `p_distance` along a Godot axis, expressed in solver space. The
// conversion is the same one positions use, because the mapping flips one sign:
// Godot +Z points along solver -Y, so the neighbour offsets are not symmetric in
// the two frames.
avbd::float3 axis_offset(int p_godot_axis, float p_distance) {
    const godot::Vector3 unit = p_godot_axis == 0 ? godot::Vector3(1, 0, 0)
                                                  : (p_godot_axis == 1 ? godot::Vector3(0, 1, 0)
                                                                       : godot::Vector3(0, 0, 1));
    return to_sim(unit * p_distance);
}

} // namespace

void AVBDSoftBody3D::set_dimensions(const Vector3i &p_dimensions) {
    dimensions = p_dimensions;
    reset_poses();
}

Vector3i AVBDSoftBody3D::get_dimensions() const {
    return dimensions;
}

void AVBDSoftBody3D::set_box_size(const Vector3 &p_size) {
    box_size = p_size;
    reset_poses();
}

Vector3 AVBDSoftBody3D::get_box_size() const {
    return box_size;
}

void AVBDSoftBody3D::set_gap(double p_gap) {
    gap = p_gap;
    reset_poses();
}

double AVBDSoftBody3D::get_gap() const {
    return gap;
}

void AVBDSoftBody3D::set_density(double p_density) {
    density = p_density;
}

double AVBDSoftBody3D::get_density() const {
    return density;
}

void AVBDSoftBody3D::set_friction(double p_friction) {
    friction = p_friction;
}

double AVBDSoftBody3D::get_friction() const {
    return friction;
}

void AVBDSoftBody3D::set_stiffness_linear(double p_stiffness) {
    stiffness_linear = p_stiffness;
}

double AVBDSoftBody3D::get_stiffness_linear() const {
    return stiffness_linear;
}

void AVBDSoftBody3D::set_stiffness_angular(double p_stiffness) {
    stiffness_angular = p_stiffness;
}

double AVBDSoftBody3D::get_stiffness_angular() const {
    return stiffness_angular;
}

int AVBDSoftBody3D::get_cell_count() const {
    return cell_counts().x * cell_counts().y * cell_counts().z;
}

// Lattice resolution with every axis clamped to at least one cell.
Vector3i AVBDSoftBody3D::cell_counts() const {
    return Vector3i(std::max(1, dimensions.x), std::max(1, dimensions.y), std::max(1, dimensions.z));
}

avbd::Rigid *AVBDSoftBody3D::_get_cell_rigid(std::size_t p_index) const {
    return p_index < rigids.size() ? rigids[p_index] : nullptr;
}

Transform3D AVBDSoftBody3D::get_cell_transform(std::size_t p_index) const {
    const avbd::Rigid *body = _get_cell_rigid(p_index);
    if (body == nullptr) {
        return Transform3D();
    }
    return sim_in_node * sim_transform(body->positionLin, body->positionAng);
}

// The solver is about to be destroyed, so `rigids` becomes dangling; `cells` keeps the
// poses for the next build.
void AVBDSoftBody3D::_clear() {
    rigids.clear();
}

void AVBDSoftBody3D::_store_state() {
    if (!cells.has_value()) {
        return;
    }
    // zip truncates to the shorter sequence, matching the min() count it replaced.
    for (auto [cell, rigid] : std::views::zip(*cells, rigids)) {
        cell.position = rigid->positionLin;
        cell.rotation = rigid->positionAng;
    }
}

// The lattice is laid out in Godot axes, centred on this node's origin, then moved
// into solver space through the node's own transform.
void AVBDSoftBody3D::_build(avbd::Solver &p_solver, const Transform3D &p_node_in_sim) {
    const Vector3i counts = cell_counts();
    const int nx = counts.x;
    const int ny = counts.y;
    const int nz = counts.z;
    const int count = nx * ny * nz;

    // Rebuilding the lattice from scratch would teleport every cell, so the pose a previous
    // build arrived at is carried over whenever the resolution is unchanged.
    const bool reuse = cells.has_value() && cells->size() == static_cast<std::size_t>(count);
    if (!reuse) {
        cells = std::vector<Cell>(count, Cell{avbd::float3{0, 0, 0}, avbd::quat::Identity()});
    }

    // Keep the MultiMesh in step with the lattice parameters.
    Ref<MultiMesh> mm = get_multimesh();
    if (mm.is_null()) {
        mm.instantiate();
        set_multimesh(mm);
    }
    if (mm->get_mesh().is_null()) {
        Ref<BoxMesh> mesh;
        mesh.instantiate();
        mesh->set_size(box_size);
        mm->set_mesh(mesh);
    }
    mm->set_transform_format(MultiMesh::TRANSFORM_3D);
    mm->set_instance_count(count);

    const avbd::float3 extents = to_sim_extents(box_size);
    const avbd::float3 half = extents * 0.5f;
    const float gap_m = static_cast<float>(gap);
    const avbd::float3 spacing = to_sim_extents(box_size + Vector3(gap_m, gap_m, gap_m));

    // Cell centres, relative to this node's origin, laid out along the Godot axes.
    const auto cell_index = [ny, nz](int x, int y, int z) { return (x * ny + y) * nz + z; };
    const auto cell_offset = [nx, ny, nz, &spacing](int x, int y, int z) {
        const float cx = static_cast<float>(x) - static_cast<float>(nx - 1) * 0.5f;
        const float cy = static_cast<float>(y) - static_cast<float>(ny - 1) * 0.5f;
        const float cz = static_cast<float>(z) - static_cast<float>(nz - 1) * 0.5f;
        return Vector3(cx * spacing.x(), cy * spacing.y(), cz * spacing.z());
    };

    rigids.assign(count, nullptr);
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                const int index = cell_index(x, y, z);

                avbd::float3 position;
                avbd::quat rotation = avbd::quat::Identity();
                if (reuse) {
                    position = (*cells)[index].position;
                    rotation = (*cells)[index].rotation;
                } else {
                    position = to_sim(p_node_in_sim.xform(cell_offset(x, y, z)));
                }

                avbd::Rigid *body = new avbd::Rigid(&p_solver, extents, static_cast<float>(density),
                        static_cast<float>(friction), position);
                body->positionAng = rotation;
                rigids[index] = body;
                (*cells)[index].position = position;
                (*cells)[index].rotation = rotation;
            }
        }
    }

    // Rigid links along every axis, plus collision suppression for the diagonal
    // neighbours that the links do not already cover (constraint-connected pairs
    // never collide, see Rigid::constrainedTo).
    const float k_linear = stiffness_from_property(stiffness_linear);
    const float k_angular = stiffness_from_property(stiffness_angular);

    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                avbd::Rigid *a = rigids[cell_index(x, y, z)];
                for (int axis = 0; axis < 3; axis++) {
                    const int cx = x + (axis == 0 ? 1 : 0);
                    const int cy = y + (axis == 1 ? 1 : 0);
                    const int cz = z + (axis == 2 ? 1 : 0);
                    if (cx >= nx || cy >= ny || cz >= nz) {
                        continue;
                    }
                    const float reach = half[SIM_AXIS[axis]];
                    new avbd::Joint(&p_solver, a, rigids[cell_index(cx, cy, cz)], axis_offset(axis, reach),
                            axis_offset(axis, -reach), k_linear, k_angular);
                }
            }
        }
    }

    for (int x = 1; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 1; z < nz; z++) {
                new avbd::IgnoreCollision(&p_solver, rigids[cell_index(x - 1, y, z - 1)], rigids[cell_index(x, y, z)]);
                new avbd::IgnoreCollision(&p_solver, rigids[cell_index(x, y, z - 1)], rigids[cell_index(x - 1, y, z)]);
            }
        }
    }
    for (int x = 0; x < nx; x++) {
        for (int y = 1; y < ny; y++) {
            for (int z = 1; z < nz; z++) {
                new avbd::IgnoreCollision(&p_solver, rigids[cell_index(x, y - 1, z - 1)], rigids[cell_index(x, y, z)]);
                new avbd::IgnoreCollision(&p_solver, rigids[cell_index(x, y, z - 1)], rigids[cell_index(x, y - 1, z)]);
            }
        }
    }
    for (int x = 1; x < nx; x++) {
        for (int y = 1; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                new avbd::IgnoreCollision(&p_solver, rigids[cell_index(x - 1, y - 1, z)], rigids[cell_index(x, y, z)]);
                new avbd::IgnoreCollision(&p_solver, rigids[cell_index(x, y - 1, z)], rigids[cell_index(x - 1, y, z)]);
            }
        }
    }
}

void AVBDSoftBody3D::_sync_visuals(const Transform3D &p_sim_in_node) {
    sim_in_node = p_sim_in_node;
    Ref<MultiMesh> mm = get_multimesh();
    if (mm.is_null()) {
        return;
    }
    const int count = std::min(mm->get_instance_count(), static_cast<int>(rigids.size()));
    for (int i = 0; i < count; i++) {
        const avbd::Rigid *body = rigids[i];
        mm->set_instance_transform(i, p_sim_in_node * sim_transform(body->positionLin, body->positionAng));
    }
}

void AVBDSoftBody3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_dimensions", "dimensions"), &AVBDSoftBody3D::set_dimensions);
    ClassDB::bind_method(D_METHOD("get_dimensions"), &AVBDSoftBody3D::get_dimensions);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "dimensions"), "set_dimensions", "get_dimensions");

    ClassDB::bind_method(D_METHOD("set_box_size", "size"), &AVBDSoftBody3D::set_box_size);
    ClassDB::bind_method(D_METHOD("get_box_size"), &AVBDSoftBody3D::get_box_size);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "box_size", PROPERTY_HINT_NONE, "suffix:m"), "set_box_size",
            "get_box_size");

    ClassDB::bind_method(D_METHOD("set_gap", "gap"), &AVBDSoftBody3D::set_gap);
    ClassDB::bind_method(D_METHOD("get_gap"), &AVBDSoftBody3D::get_gap);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "gap", PROPERTY_HINT_RANGE, "0,10,0.001,or_greater,suffix:m"), "set_gap",
            "get_gap");

    ClassDB::bind_method(D_METHOD("set_density", "density"), &AVBDSoftBody3D::set_density);
    ClassDB::bind_method(D_METHOD("get_density"), &AVBDSoftBody3D::get_density);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "density", PROPERTY_HINT_RANGE, "0,10000,0.01,or_greater"), "set_density",
            "get_density");

    ClassDB::bind_method(D_METHOD("set_friction", "friction"), &AVBDSoftBody3D::set_friction);
    ClassDB::bind_method(D_METHOD("get_friction"), &AVBDSoftBody3D::get_friction);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "friction", PROPERTY_HINT_RANGE, "0,1,0.01,or_greater"), "set_friction",
            "get_friction");

    ClassDB::bind_method(D_METHOD("set_stiffness_linear", "stiffness"), &AVBDSoftBody3D::set_stiffness_linear);
    ClassDB::bind_method(D_METHOD("get_stiffness_linear"), &AVBDSoftBody3D::get_stiffness_linear);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "stiffness_linear", PROPERTY_HINT_RANGE, "-1,1000000,1,or_greater"),
            "set_stiffness_linear", "get_stiffness_linear");

    ClassDB::bind_method(D_METHOD("set_stiffness_angular", "stiffness"), &AVBDSoftBody3D::set_stiffness_angular);
    ClassDB::bind_method(D_METHOD("get_stiffness_angular"), &AVBDSoftBody3D::get_stiffness_angular);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "stiffness_angular", PROPERTY_HINT_RANGE, "-1,100000,1,or_greater"),
            "set_stiffness_angular", "get_stiffness_angular");

    ClassDB::bind_method(D_METHOD("get_cell_count"), &AVBDSoftBody3D::get_cell_count);
    ClassDB::bind_method(D_METHOD("get_cell_transform", "index"), &AVBDSoftBody3D::get_cell_transform);
}

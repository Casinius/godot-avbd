/*
 * AVBDSoftBody3D - a lattice of jointed boxes that behaves like a soft body.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * This is the formulation the AVBD authors use for soft bodies: a W x D x H grid of
 * rigid boxes wired together with stiff joints, with diagonal neighbours whose
 * collisions are suppressed. The boxes are real solver bodies, so the lattice
 * collides with the rest of the scene (and with itself along its outer shell). The
 * lattice is drawn through a single MultiMesh, so a 4x4x4 blob costs one draw call
 * and one node.
 *
 * The lattice is centred on this node's origin: the node's transform places it.
 */

#ifndef AVBD_SOFT_BODY3D_HPP
#define AVBD_SOFT_BODY3D_HPP

#include <optional>
#include <cstddef>
#include <vector>

#include <godot_cpp/classes/multi_mesh_instance3d.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/vector3.hpp>
#include <godot_cpp/variant/vector3i.hpp>

#include "avbd/solver.h"

namespace godot {

class AVBDSoftBody3D : public MultiMeshInstance3D {
    GDCLASS(AVBDSoftBody3D, MultiMeshInstance3D)

    Vector3i dimensions = Vector3i(4, 4, 4);
    Vector3 box_size = Vector3(0.25, 0.25, 0.25);
    double gap = 0.0;
    double density = 1.0;
    double friction = 0.5;
    double stiffness_linear = 1000.0;
    double stiffness_angular = 250.0;

    // Latest simulated pose of each cell. A rebuild recreates the solver, so the poses are
    // kept here to put the lattice back where it was rather than resetting it to its spawn
    // layout. Absent until the lattice has been built once.
    struct Cell {
        avbd::float3 position;
        avbd::quat rotation;
    };

    std::vector<avbd::Rigid *> rigids;
    std::optional<std::vector<Cell>> cells;
    // Solver space -> this node's local space, refreshed every simulated frame.
    Transform3D sim_in_node;

protected:
    static void _bind_methods();

    // Drop the cached cell poses, so the next build lays the lattice out from scratch.
    void reset_poses() { cells.reset(); }

public:
    AVBDSoftBody3D() = default;
    ~AVBDSoftBody3D() override = default;

    // Lattice resolution along each Godot axis.
    void set_dimensions(const Vector3i &p_dimensions);
    Vector3i get_dimensions() const;

    // Size of each box in the lattice (full extents, Godot axes).
    void set_box_size(const Vector3 &p_size);
    Vector3 get_box_size() const;

    // Extra spacing between neighbouring boxes; 0 makes them touch.
    void set_gap(double p_gap);
    double get_gap() const;

    void set_density(double p_density);
    double get_density() const;
    void set_friction(double p_friction);
    double get_friction() const;

    // Joint stiffness between neighbouring boxes: negative is a rigid link.
    void set_stiffness_linear(double p_stiffness);
    double get_stiffness_linear() const;
    void set_stiffness_angular(double p_stiffness);
    double get_stiffness_angular() const;

    int get_cell_count() const;
    // Lattice resolution, with every axis clamped to at least one cell.
    Vector3i cell_counts() const;

    // Pose of one lattice cell, in this node's local space. Cells are indexed
    // x-major, then y, then z (x * ny * nz + y * nz + z).
    Transform3D get_cell_transform(std::size_t p_index) const;

    // --- simulation interface, used by AVBDWorld3D ---------------------------
    void _build(avbd::Solver &p_solver, const Transform3D &p_node_in_sim);
    void _store_state();
    void _clear();
    void _sync_visuals(const Transform3D &p_sim_in_node);
    avbd::Rigid *_get_cell_rigid(std::size_t p_index) const;
};

} // namespace godot

#endif // AVBD_SOFT_BODY3D_HPP

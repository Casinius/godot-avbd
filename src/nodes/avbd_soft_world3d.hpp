/*
 * AVBDSoftWorld3D - drives AVBDSoftBody3D nodes through a solver of their own.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * This is the minimal surviving piece of the old AVBDWorld3D node layer: soft bodies
 * have no PhysicsServer3D representation this round, so a node has to own their solver.
 * Rigid bodies and joints do NOT go here - a scene with those uses the AVBD physics
 * server (project setting physics/3d/physics_engine = "AVBD") and standard Godot nodes.
 *
 * Responsibilities, deliberately narrow:
 *   - collect the AVBDSoftBody3D nodes in this subtree (re-scan on tree changes)
 *   - collect child StaticBody3D "ground" nodes (one box shape each) as static solver
 *     bodies, so a lattice has something to rest on: soft cells live in this node's
 *     solver and cannot collide with bodies the physics server owns
 *   - step their solver in _physics_process (following the game: a paused SceneTree
 *     pauses this node through the engine's own pause handling)
 *   - write results back to the soft body visuals
 * Properties: gravity and threads only. Everything else is the soft body's own.
 */

#ifndef AVBD_SOFT_WORLD3D_HPP
#define AVBD_SOFT_WORLD3D_HPP

#include <vector>

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/box_shape3d.hpp>
#include <godot_cpp/classes/collision_shape3d.hpp>
#include <godot_cpp/classes/static_body3d.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

namespace godot {

class AVBDSoftBody3D;

class AVBDSoftWorld3D : public Node3D {
    GDCLASS(AVBDSoftWorld3D, Node3D)

    double gravity = 9.8;
    int threads = 0;

    avbd::Solver solver;
    std::vector<AVBDSoftBody3D *> soft_bodies;

    // Static ground boxes, mirrored into the solver so lattices can rest on them.
    struct GroundEntry {
        StaticBody3D *node = nullptr;
        avbd::Rigid *rigid = nullptr;
        Vector3 extents; // full extents, Godot axes
    };
    std::vector<GroundEntry> grounds;

    // Scratch buffer for the per-tick subtree scan.
    std::vector<AVBDSoftBody3D *> scan_soft_bodies;
    std::vector<GroundEntry> scan_grounds;

protected:
    static void _bind_methods();
    void _notification(int p_what);

    void _collect(Node *p_node);
    bool _scan_changed();
    static Shape3D *_first_shape(StaticBody3D *p_body);
public:
    AVBDSoftWorld3D() = default;
    ~AVBDSoftWorld3D() override;

    // Gravity in m/s^2, applied along -Y (Godot's down direction).
    void set_gravity(double p_gravity);
    double get_gravity() const;

    // Worker threads for the solver's per-body phases. 0 = one per hardware thread,
    // 1 = run inline. The result does not depend on this value.
    void set_threads(int p_threads);
    int get_threads() const;

    // Number of lattice cells currently simulated (reporting only).
    int get_body_count() const;

    // Force a re-scan and rebuild without waiting for the next tree change.
    void rebuild();

    // Step the solver for one physics tick.
    void _physics_process(double p_delta) override;

    // Simulation coordinate system: this node's transform (as before).
    void _rebuild();
    void _release();
    void _ready() override;
};

} // namespace godot

#endif // AVBD_SOFT_WORLD3D_HPP

/* AVBDSoftWorld3D drives AVBDSoftBody3D nodes in their own solver (no server representation);
 * rigid bodies/joints use the AVBD physics server. Core ported from savant117/avbd-demo3d (MIT). */

#ifndef AVBD_SOFT_WORLD3D_HPP
#define AVBD_SOFT_WORLD3D_HPP

#include <vector>

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/box_shape3d.hpp>
#include <godot_cpp/classes/collision_shape3d.hpp>
#include <godot_cpp/classes/static_body3d.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "avbd/solver.h"

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

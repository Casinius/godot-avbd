/*
 * AVBDConstraint3D - shared base for every AVBD force that connects two bodies.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#ifndef AVBD_CONSTRAINT3D_HPP
#define AVBD_CONSTRAINT3D_HPP

#include <unordered_map>

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/variant/node_path.hpp>
#include <godot_cpp/variant/transform3d.hpp>

#include "avbd/solver.h"

namespace godot {

class AVBDRigidBody3D;

// Bodies taking part in one simulation, looked up by their scene node.
using AVBDBodyMap = std::unordered_map<Node3D *, avbd::Rigid *>;

class AVBDConstraint3D : public Node3D {
    GDCLASS(AVBDConstraint3D, Node3D)

    NodePath node_a;
    NodePath node_b;

protected:
    static void _bind_methods();

    // The solver-side force this node created, or null when not simulated.
    avbd::Force *force = nullptr;
    // Set once per tick by the world: whether `force` is still in the solver.
    bool force_alive = false;

    // Called when the solver removed this node's force (a broken joint is deleted by
    // the solver, which also stops it from suppressing collisions between its bodies).
    virtual void _on_force_removed() {}

public:
    AVBDConstraint3D() = default;
    ~AVBDConstraint3D() override = default;

    // Body A. Leave empty for joints that pin body B to a fixed point in space.
    void set_node_a(const NodePath &p_path);
    NodePath get_node_a() const;

    // Body B. Always required.
    void set_node_b(const NodePath &p_path);
    NodePath get_node_b() const;

    AVBDRigidBody3D *get_body_a() const;
    AVBDRigidBody3D *get_body_b() const;

    // True once the owning world has created this constraint in the solver. False
    // means the endpoints could not be resolved, or the node is outside any world.
    bool is_simulated() const;

    // Create the solver-side force, or null if the endpoints cannot be resolved.
    // `p_world_global` is the transform of the owning AVBDWorld3D, used to interpret
    // world-space anchors. Subclasses override this; the base returns null because
    // godot-cpp needs every registered class to be instantiable.
    virtual avbd::Force *create_force(avbd::Solver &, const AVBDBodyMap &, const Transform3D &) const {
        return nullptr;
    }

    // --- simulation interface, used by AVBDWorld3D ---------------------------
    // The world keeps the force it got from create_force() and hands it back here, so
    // a constraint node can report solver-side state (e.g. whether it has broken).
    void _bind_force(avbd::Force *p_force);
    bool _is_simulated() const { return force != nullptr; }
    // Whether the constraint is permanently inactive (a broken joint): the world then
    // creates no solver force for it, which also lets its bodies collide again.
    virtual bool _is_inactive() const { return false; }
    // Per-tick force liveness check, driven by AVBDWorld3D.
    void _mark_force_alive(bool p_alive) { force_alive = p_alive; }
    void _resolve_force_removal();
};

} // namespace godot

#endif // AVBD_CONSTRAINT3D_HPP

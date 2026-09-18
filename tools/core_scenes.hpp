/*
 * Standard AVBD test scenes, ported from the authors' reference implementation
 * (github.com/savant117/avbd-demo3d, MIT, Copyright (c) 2026 Chris Giles).
 *
 * Units: metres, kilograms, seconds. Gravity acts along -Z; +Z is up.
 */

#ifndef AVBD_CORE_SCENES_HPP
#define AVBD_CORE_SCENES_HPP

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {

struct CoreScene {
    const char *name;
    void (*build)(Solver *solver);
};

// Named scenes, in the reference implementation's order.
extern const CoreScene coreScenes[];
extern const int coreSceneCount;

} // namespace avbd

#endif // AVBD_CORE_SCENES_HPP

/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#pragma once

#include <cmath>
#include <cstdint>

namespace avbd {

// Converts a stream of irregular engine steps into whole simulation ticks.
//
// The solver's stability belongs to its own timestep, not to whatever rate the engine happens
// to call at, so the engine step goes into an accumulator and only whole ticks of logical time
// come out. Two rules keep that honest:
//
//  - The only time source is the step the engine passes in. Reading a wall clock here would
//    make the tick sequence depend on host timing, and the solver's determinism contract (same
//    result for any thread count, in any process) would be about elapsed real seconds.
//  - Debt is dropped, never carried. A slow frame buys `max_ticks` ticks and the surplus logical
//    time is discarded, which dilates time instead of enlarging the internal dt. A solver step
//    always covers exactly `tick_dt`, so a catch-up never simulates a longer step than the
//    steady state does.
//
// The leftover under one tick is reported as `r_alpha`: the fraction of the next tick that the
// caller has already accumulated, which is what a renderer interpolates with.
struct FixedClock {
    double tick_dt = 1.0 / 120.0;
    int max_ticks = 4;

    double acc = 0.0;          // logical time in hand but not yet simulated
    uint64_t dropped_ticks = 0; // whole ticks discarded by the catch-up cap

    // Consumes `p_step` and returns the number of ticks to run now (0 is a valid answer: an
    // engine step shorter than one tick runs nothing and leaves the time in `acc`).
    int consume(double p_step, double &r_alpha) {
        if (!(tick_dt > 0.0)) { // an unset rate steps no time at all; keep the bank empty
            acc = 0.0;
            r_alpha = 0.0;
            return 0;
        }
        acc += p_step;
        // The cap is applied on the double, before the conversion: a pathological step (a
        // debugger left open, a level load) can bank hundreds of ticks and the count of dropped
        // ones has to stay meaningful.
        const double affordable = acc / tick_dt;
        if (affordable > static_cast<double>(max_ticks)) {
            acc = std::fmod(acc, tick_dt); // the debt is dropped, so nothing carries it forward
            dropped_ticks += static_cast<uint64_t>(affordable) - static_cast<uint64_t>(max_ticks);
            r_alpha = acc / tick_dt;
            return max_ticks;
        }
        const int ticks = affordable > 0.0 ? static_cast<int>(affordable) : 0;
        acc -= static_cast<double>(ticks) * tick_dt;
        r_alpha = acc / tick_dt;
        return ticks;
    }

    // Called when stepping restarts (the server or a space was inactive): time that passed
    // while nothing was simulated is not a debt the next frame has to pay.
    void reset() {
        acc = 0.0;
    }
};

} // namespace avbd

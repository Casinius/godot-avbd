/*
 * Shared broad-phase pair encoding.
 *
 * Both the hierarchical hash and the Morton LBVH emit the same packed candidate key, so
 * the sort/dedupe/filter downstream is one implementation rather than two. The key is
 * (smaller body index << 32) | larger body index: lexicographic order over std::pair is
 * plain integer order, which is what lets the pair list use a radix sort.
 */

#pragma once

#include <cstdint>

namespace avbd {

[[nodiscard]] inline constexpr uint64_t packPair(int p_a, int p_b) noexcept
{
    const uint64_t lo = static_cast<uint64_t>(static_cast<uint32_t>(p_a < p_b ? p_a : p_b));
    const uint64_t hi = static_cast<uint64_t>(static_cast<uint32_t>(p_a < p_b ? p_b : p_a));
    return (lo << 32) | hi;
}

[[nodiscard]] inline constexpr int pairFirst(uint64_t p_key) noexcept
{
    return static_cast<int>(p_key >> 32);
}

[[nodiscard]] inline constexpr int pairSecond(uint64_t p_key) noexcept
{
    return static_cast<int>(p_key & 0xFFFFFFFFu);
}

} // namespace avbd

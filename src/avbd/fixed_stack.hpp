/*
 * FixedStack - a bounded LIFO for hot traversals that must not allocate.
 *
 * The BVH descents (broad-phase self-query, pick, AABB queries) walk trees of known
 * maximum depth, where a raw `int stack[64]` plus a hand-kept index works until the
 * day the depth bound changes and the overflow silently corrupts the caller's frame.
 * This is that array with the two failure modes turned into assertions (compiled out
 * under NDEBUG, so the release path is the same instructions as the raw array) and
 * the index bookkeeping in one place.
 */

#ifndef AVBD_FIXED_STACK_HPP
#define AVBD_FIXED_STACK_HPP

#include <array>
#include <cassert>
#include <cstddef>

namespace avbd {

template <typename T, std::size_t kCapacity>
class FixedStack {
public:
    void push(const T &value) {
        assert(size_ < kCapacity && "FixedStack overflow: raise kCapacity");
        data_[size_] = value;
        ++size_;
    }

    [[nodiscard]] T pop() {
        assert(size_ > 0 && "FixedStack underflow: pop on empty");
        --size_;
        return data_[size_];
    }

    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
    [[nodiscard]] std::size_t size() const noexcept { return size_; }

private:
    std::array<T, kCapacity> data_{};
    std::size_t size_ = 0;
};

} // namespace avbd

#endif // AVBD_FIXED_STACK_HPP

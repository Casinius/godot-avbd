// radix_sort.hpp —— 与 std::sort 接口一致的基数排序后端
#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <type_traits>
#include <vector>

namespace radix_std {

namespace detail {

// 仅对这些整数类型启用基数排序
template <typename T> struct is_radix_candidate : std::false_type {};
template <> struct is_radix_candidate<std::int32_t>  : std::true_type {};
template <> struct is_radix_candidate<std::uint32_t> : std::true_type {};
template <> struct is_radix_candidate<std::int64_t>  : std::true_type {};
template <> struct is_radix_candidate<std::uint64_t> : std::true_type {};

template <typename T>
bool radix_sort_bytes(T* data, std::size_t n) {
    using U = std::make_unsigned_t<T>;
    constexpr bool is_signed = std::is_signed_v<T>;
    constexpr unsigned bits = sizeof(U) * 8;

    // 把元素解释成"用于排序的无符号位模式"：有符号数翻转符号位，使负数排在正数之前。
    // 使用 memcpy 而非 reinterpret_cast，避免 strict-aliasing 问题，也不依赖字节序。
    auto key_of = [](const T& v) -> U {
        U u;
        std::memcpy(&u, &v, sizeof(U));
        if constexpr (is_signed) u ^= (U(1) << (bits - 1));
        return u;
    };

    // 分配辅助缓冲；失败则返回 false 让调用者回退
    std::vector<T> buf;
    try { buf.resize(n); }
    catch (...) { return false; }

    T* src = data;
    T* dst = buf.data();
    bool swapped = false;
    std::array<std::size_t, 256> count{};

    // LSD：每 8 位为一趟，从最低位开始
    for (unsigned shift = 0; shift < bits; shift += 8) {
        count.fill(0);
        for (std::size_t i = 0; i < n; ++i)
            ++count[(key_of(src[i]) >> shift) & 0xFFu];

        std::size_t total = 0;
        for (auto& c : count) { auto old = c; c = total; total += old; }

        for (std::size_t i = 0; i < n; ++i) {
            auto b = (key_of(src[i]) >> shift) & 0xFFu;
            dst[count[b]++] = src[i];
        }
        std::swap(src, dst);
        swapped = !swapped;
    }
    if (swapped) std::copy(src, src + n, data);
    return true;
}

} // namespace detail

// 三参数版本：仅当 Compare 正好是 std::less<T> 且类型是 32/64 位整数时尝试基数排序
template <typename RandomIt, typename Compare>
void sort(RandomIt first, RandomIt last, Compare comp) {
    using T = typename std::iterator_traits<RandomIt>::value_type;

    if constexpr (std::is_same_v<Compare, std::less<T>> &&
                  detail::is_radix_candidate<T>::value) {
        const auto n = static_cast<std::size_t>(last - first);
        // 小规模时比较排序常数更小，直接回退
        if (n >= 2 && n >= 64 && detail::radix_sort_bytes(&*first, n)) return;
    }
    std::sort(first, last, comp);
}

// 两参数版本：与 std::sort 一致的默认比较器
template <typename RandomIt>
void sort(RandomIt first, RandomIt last) {
    using T = typename std::iterator_traits<RandomIt>::value_type;
    radix_std::sort(first, last, std::less<T>{});
}

} // namespace radix_std
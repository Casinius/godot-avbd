/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 *
 * Port note: an earlier revision of this file contained `using namespace std;`, which
 * leaked every standard-library name into any translation unit that included it. The
 * handful of places that relied on it now name `std::` explicitly.
 */

#pragma once

#include <cmath>
#include <cstddef>

namespace avbd {

// -----------------------------------------------------------------------------
// Math types
//
// Plain aggregates: no constructors, no invariants to break, so they stay trivially
// copyable and cheap to pass by value. Component access goes through `operator[]`, which
// spells the field out rather than aliasing the struct's address as an array - that trick
// is undefined behaviour, and it was the only reason this header needed the C-style casts
// it used to contain. With a constant index the comparison folds away entirely.
// -----------------------------------------------------------------------------

struct float2
{
    float x, y;

    float &operator[](std::size_t i) noexcept { return i == 0 ? x : y; }
    const float &operator[](std::size_t i) const noexcept { return i == 0 ? x : y; }
};

struct float3
{
    float x, y, z;

    float &operator[](std::size_t i) noexcept { return i == 0 ? x : (i == 1 ? y : z); }
    const float &operator[](std::size_t i) const noexcept { return i == 0 ? x : (i == 1 ? y : z); }
};

struct quat
{
    float x, y, z, w;

    float &operator[](std::size_t i) noexcept { return i == 0 ? x : (i == 1 ? y : (i == 2 ? z : w)); }
    const float &operator[](std::size_t i) const noexcept { return i == 0 ? x : (i == 1 ? y : (i == 2 ? z : w)); }
};

struct float2x2
{
    float2 row[2];

    float2 &operator[](std::size_t i) noexcept { return row[i]; }
    const float2 &operator[](std::size_t i) const noexcept { return row[i]; }

    [[nodiscard]] float2 col(std::size_t i) const noexcept { return float2{row[0][i], row[1][i]}; }
};

struct float3x3
{
    float3 row[3];

    float3 &operator[](std::size_t i) noexcept { return row[i]; }
    const float3 &operator[](std::size_t i) const noexcept { return row[i]; }

    [[nodiscard]] float3 col(std::size_t i) const noexcept { return float3{row[0][i], row[1][i], row[2][i]}; }
};

// Declared up front: the quaternion operators below are defined in terms of these.
[[nodiscard]] quat normalize(quat q) noexcept;
[[nodiscard]] quat inverse(quat q) noexcept;
[[nodiscard]] quat conjugate(quat q) noexcept;
[[nodiscard]] quat operator+(quat a, float3 b) noexcept;

// -----------------------------------------------------------------------------
// float2 operators
// -----------------------------------------------------------------------------

inline float2 &operator+=(float2 &a, float2 b) noexcept
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

inline float2 &operator-=(float2 &a, float2 b) noexcept
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

[[nodiscard]] inline float2 operator-(float2 v) noexcept
{
    return {-v.x, -v.y};
}

[[nodiscard]] inline float2 operator+(float2 a, float2 b) noexcept
{
    return {a.x + b.x, a.y + b.y};
}

[[nodiscard]] inline float2 operator-(float2 a, float2 b) noexcept
{
    return {a.x - b.x, a.y - b.y};
}

[[nodiscard]] inline float2 operator*(float2 a, float b) noexcept
{
    return {a.x * b, a.y * b};
}

[[nodiscard]] inline float2 operator/(float2 a, float b) noexcept
{
    return {a.x / b, a.y / b};
}

// -----------------------------------------------------------------------------
// float3 operators
// -----------------------------------------------------------------------------

inline float3 &operator+=(float3 &a, float3 b) noexcept
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline float3 &operator-=(float3 &a, float3 b) noexcept
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

[[nodiscard]] inline float3 operator-(float3 v) noexcept
{
    return {-v.x, -v.y, -v.z};
}

[[nodiscard]] inline float3 operator+(float3 a, float3 b) noexcept
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

[[nodiscard]] inline float3 operator-(float3 a, float3 b) noexcept
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

[[nodiscard]] inline float3 operator*(float3 a, float b) noexcept
{
    return {a.x * b, a.y * b, a.z * b};
}

[[nodiscard]] inline float3 operator/(float3 a, float b) noexcept
{
    return {a.x / b, a.y / b, a.z / b};
}

// -----------------------------------------------------------------------------
// float2x2 operators
// -----------------------------------------------------------------------------

[[nodiscard]] inline float2x2 operator+(float2x2 a, float2x2 b) noexcept
{
    return {a[0] + b[0], a[1] + b[1]};
}

[[nodiscard]] inline float2x2 operator-(float2x2 a, float2x2 b) noexcept
{
    return {a[0] - b[0], a[1] - b[1]};
}

[[nodiscard]] inline float2x2 operator*(float2x2 a, float b) noexcept
{
    return {a[0] * b, a[1] * b};
}

[[nodiscard]] inline float2x2 operator/(float2x2 a, float b) noexcept
{
    return {a[0] / b, a[1] / b};
}

// -----------------------------------------------------------------------------
// float3x3 operators
// -----------------------------------------------------------------------------

inline float3x3 &operator+=(float3x3 &a, float3x3 b) noexcept
{
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    return a;
}

[[nodiscard]] inline float3x3 operator+(float3x3 a, float3x3 b) noexcept
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

[[nodiscard]] inline float3x3 operator-(float3x3 a, float3x3 b) noexcept
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

[[nodiscard]] inline float3x3 operator*(float3x3 a, float b) noexcept
{
    return {a[0] * b, a[1] * b, a[2] * b};
}

[[nodiscard]] inline float3x3 operator/(float3x3 a, float b) noexcept
{
    return {a[0] / b, a[1] / b, a[2] / b};
}

[[nodiscard]] inline float3x3 operator-(float3x3 a) noexcept
{
    return {-a.row[0], -a.row[1], -a.row[2]};
}

// -----------------------------------------------------------------------------
// quat operators
// -----------------------------------------------------------------------------

[[nodiscard]] inline quat operator*(quat a, float b) noexcept
{
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}

[[nodiscard]] inline quat operator/(quat a, float b) noexcept
{
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}

[[nodiscard]] inline quat operator*(quat a, quat b) noexcept
{
    return {
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

[[nodiscard]] inline quat operator+(quat a, quat b) noexcept
{
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

// Vector part of the rotation that takes `b` to `a`, scaled by two: the small-angle
// rotation difference the angular constraint is built on.
[[nodiscard]] inline float3 operator-(quat a, quat b) noexcept
{
    const quat d = a * inverse(b);
    return float3{d.x, d.y, d.z} * 2.0f;
}

[[nodiscard]] inline quat operator+(quat a, float3 b) noexcept
{
    return normalize(a + quat{b.x, b.y, b.z, 0} * a * 0.5f);
}

// -----------------------------------------------------------------------------
// Scalar and vector functions
//
// `abs(float)` and friends are declared here on purpose: without them, unqualified
// `abs(x)` on a float would not resolve, since the vector overloads below hide the
// C library's. They are exact replacements for std::fabs/std::abs on a float.
// -----------------------------------------------------------------------------

[[nodiscard]] inline float abs(float x) noexcept
{
    return std::fabs(x);
}

[[nodiscard]] inline float2 abs(float2 v) noexcept
{
    return {std::fabs(v.x), std::fabs(v.y)};
}

[[nodiscard]] inline float3 abs(float3 v) noexcept
{
    return {std::fabs(v.x), std::fabs(v.y), std::fabs(v.z)};
}

[[nodiscard]] inline float2x2 abs(float2x2 a) noexcept
{
    return {abs(a[0]), abs(a[1])};
}

[[nodiscard]] inline float rad(float deg) noexcept
{
    return deg * 0.01745329251994329577f;
}

[[nodiscard]] inline float sign(float x) noexcept
{
    return x < 0 ? -1.0f : x > 0 ? 1.0f
                                 : 0.0f;
}

[[nodiscard]] inline float min(float a, float b) noexcept
{
    return a < b ? a : b;
}

[[nodiscard]] inline float max(float a, float b) noexcept
{
    return a > b ? a : b;
}

[[nodiscard]] inline float3 min(float3 a, float b) noexcept
{
    return {min(a.x, b), min(a.y, b), min(a.z, b)};
}

[[nodiscard]] inline float clamp(float x, float a, float b) noexcept
{
    return max(a, min(b, x));
}

[[nodiscard]] inline float3 clamp(float3 v, float a, float b) noexcept
{
    return {clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b)};
}

[[nodiscard]] inline float dot(float2 a, float2 b) noexcept
{
    return a.x * b.x + a.y * b.y;
}

[[nodiscard]] inline float dot(float3 a, float3 b) noexcept
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] inline float lengthSq(float2 v) noexcept
{
    return dot(v, v);
}

[[nodiscard]] inline float length(float2 v) noexcept
{
    return std::sqrt(lengthSq(v));
}

[[nodiscard]] inline float lengthSq(float3 v) noexcept
{
    return dot(v, v);
}

[[nodiscard]] inline float length(float3 v) noexcept
{
    return std::sqrt(lengthSq(v));
}

[[nodiscard]] inline float lengthSq(quat q) noexcept
{
    return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
}

[[nodiscard]] inline float length(quat q) noexcept
{
    return std::sqrt(lengthSq(q));
}

[[nodiscard]] inline float3 normalize(float3 v) noexcept
{
    return v / length(v);
}

[[nodiscard]] inline quat normalize(quat q) noexcept
{
    return q / length(q);
}

[[nodiscard]] inline float cross(float2 a, float2 b) noexcept
{
    return a.x * b.y - a.y * b.x;
}

[[nodiscard]] inline float3 cross(float3 a, float3 b) noexcept
{
    return float3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

[[nodiscard]] inline float3x3 skew(float3 r) noexcept
{
    return float3x3{{
        0, -r.z, r.y,
        r.z, 0, -r.x,
        -r.y, r.x, 0}};
}

[[nodiscard]] inline float2x2 outer(float2 a, float2 b) noexcept
{
    return {b * a.x, b * a.y};
}

[[nodiscard]] inline float3x3 outer(float3 a, float3 b) noexcept
{
    return {b * a.x, b * a.y, b * a.z};
}

[[nodiscard]] inline float2x2 transpose(float2x2 a) noexcept
{
    return {float2{a[0][0], a[1][0]}, float2{a[0][1], a[1][1]}};
}

[[nodiscard]] inline float3x3 transpose(float3x3 a) noexcept
{
    return {float3{a[0][0], a[1][0], a[2][0]}, float3{a[0][1], a[1][1], a[2][1]}, float3{a[0][2], a[1][2], a[2][2]}};
}

[[nodiscard]] inline float3x3 diagonal(float m00, float m11, float m22) noexcept
{
    return float3x3{{
        m00, 0, 0,
        0, m11, 0,
        0, 0, m22}};
}

[[nodiscard]] inline quat conjugate(quat q) noexcept
{
    return {-q.x, -q.y, -q.z, q.w};
}

[[nodiscard]] inline quat inverse(quat q) noexcept
{
    return conjugate(q) / lengthSq(q);
}

// Rotate `v` by `angle`: the usual u v u* form, written with two cross products.
[[nodiscard]] inline float3 rotate(quat angle, float3 v) noexcept
{
    const float3 u = {angle.x, angle.y, angle.z};
    const float3 t = cross(u, v) * 2.0f;
    return v + t * angle.w + cross(u, t);
}

[[nodiscard]] inline float3 transform(float3 qLin, quat qAng, float3 v) noexcept
{
    return rotate(qAng, v) + qLin;
}

[[nodiscard]] inline float3x3 rotation(quat q) noexcept
{
    const float x = q.x, y = q.y, z = q.z, w = q.w;

    const float xx = x * x, yy = y * y, zz = z * z;
    const float xy = x * y, xz = x * z, yz = y * z;
    const float wx = w * x, wy = w * y, wz = w * z;

    // Row-major 3x3
    return float3x3{{
        1.0f - 2.0f * (yy + zz), 2.0f * (xy + wz), 2.0f * (xz - wy),
        2.0f * (xy - wz), 1.0f - 2.0f * (xx + zz), 2.0f * (yz + wx),
        2.0f * (xz + wy), 2.0f * (yz - wx), 1.0f - 2.0f * (xx + yy)}};
}

[[nodiscard]] inline float3x3 orthonormal(float3 normal) noexcept
{
    float3 t1 = std::fabs(normal.x) > std::fabs(normal.z) ? float3{-normal.y, normal.x, 0} : float3{0, -normal.z, normal.y};
    t1 = normalize(t1);
    const float3 t2 = cross(normal, t1);
    return float3x3{{normal, t1, t2}};
}

[[nodiscard]] inline float3x3 diagonalize(float3x3 m) noexcept
{
    return diagonal(length(m.col(0)), length(m.col(1)), length(m.col(2)));
}

[[nodiscard]] inline float2x2 operator*(float2x2 a, float2 b) noexcept
{
    return {dot(a[0], b), dot(a[1], b)};
}

[[nodiscard]] inline float3 operator*(float3x3 a, float3 b) noexcept
{
    return {dot(a[0], b), dot(a[1], b), dot(a[2], b)};
}

[[nodiscard]] inline float2x2 operator*(float2x2 a, float2x2 b) noexcept
{
    return {
        float2{dot(a.row[0], b.col(0)), dot(a.row[0], b.col(1))},
        float2{dot(a.row[1], b.col(0)), dot(a.row[1], b.col(1))}};
}

[[nodiscard]] inline float3x3 operator*(float3x3 a, float3x3 b) noexcept
{
    return {
        float3{dot(a.row[0], b.col(0)), dot(a.row[0], b.col(1)), dot(a.row[0], b.col(2))},
        float3{dot(a.row[1], b.col(0)), dot(a.row[1], b.col(1)), dot(a.row[1], b.col(2))},
        float3{dot(a.row[2], b.col(0)), dot(a.row[2], b.col(1)), dot(a.row[2], b.col(2))}};
}

// Solve the symmetric 6x6 system [aLin aCross^T; aCross aAng] x = [bLin; bAng] by LDL^T
// and write the result into `xLin`/`xAng`. The matrix is stored as its lower triangle:
// only the lower half of each 3x3 block is read, and the entries that are never touched
// (the strictly upper part of aLin, aAng and aCross) are ignored.
//
// This is the innermost loop of the solver: it runs once per body per iteration, so it is
// deliberately written out flat rather than expressed with the matrix types.
inline void solve(float3x3 aLin, float3x3 aAng, float3x3 aCross, float3 bLin, float3 bAng, float3 &xLin, float3 &xAng) noexcept
{
    // Extract elements from lower triangle storage
    float A11 = aLin[0][0];
    float A21 = aLin[1][0], A22 = aLin[1][1];
    float A31 = aLin[2][0], A32 = aLin[2][1], A33 = aLin[2][2];
    float A41 = aCross[0][0], A42 = aCross[0][1], A43 = aCross[0][2], A44 = aAng[0][0];
    float A51 = aCross[1][0], A52 = aCross[1][1], A53 = aCross[1][2], A54 = aAng[1][0], A55 = aAng[1][1];
    float A61 = aCross[2][0], A62 = aCross[2][1], A63 = aCross[2][2], A64 = aAng[2][0], A65 = aAng[2][1], A66 = aAng[2][2];

    // Step 1: LDL^T decomposition
    float L21 = A21 / A11;
    float L31 = A31 / A11;
    float L41 = A41 / A11;
    float L51 = A51 / A11;
    float L61 = A61 / A11;

    float D1 = A11;

    float D2 = A22 - L21 * L21 * D1;

    float L32 = (A32 - L21 * L31 * D1) / D2;
    float L42 = (A42 - L21 * L41 * D1) / D2;
    float L52 = (A52 - L21 * L51 * D1) / D2;
    float L62 = (A62 - L21 * L61 * D1) / D2;

    float D3 = A33 - (L31 * L31 * D1 + L32 * L32 * D2);

    float L43 = (A43 - L31 * L41 * D1 - L32 * L42 * D2) / D3;
    float L53 = (A53 - L31 * L51 * D1 - L32 * L52 * D2) / D3;
    float L63 = (A63 - L31 * L61 * D1 - L32 * L62 * D2) / D3;

    float D4 = A44 - (L41 * L41 * D1 + L42 * L42 * D2 + L43 * L43 * D3);

    float L54 = (A54 - L41 * L51 * D1 - L42 * L52 * D2 - L43 * L53 * D3) / D4;
    float L64 = (A64 - L41 * L61 * D1 - L42 * L62 * D2 - L43 * L63 * D3) / D4;

    float D5 = A55 - (L51 * L51 * D1 + L52 * L52 * D2 + L53 * L53 * D3 + L54 * L54 * D4);

    float L65 = (A65 - L51 * L61 * D1 - L52 * L62 * D2 - L53 * L63 * D3 - L54 * L64 * D4) / D5;

    float D6 = A66 - (L61 * L61 * D1 + L62 * L62 * D2 + L63 * L63 * D3 + L64 * L64 * D4 + L65 * L65 * D5);

    // Step 2: Forward substitution: Solve Ly = b
    float y1 = bLin[0];
    float y2 = bLin[1] - L21 * y1;
    float y3 = bLin[2] - L31 * y1 - L32 * y2;
    float y4 = bAng[0] - L41 * y1 - L42 * y2 - L43 * y3;
    float y5 = bAng[1] - L51 * y1 - L52 * y2 - L53 * y3 - L54 * y4;
    float y6 = bAng[2] - L61 * y1 - L62 * y2 - L63 * y3 - L64 * y4 - L65 * y5;

    // Step 3: Diagonal solve: Solve Dz = y
    float z1 = y1 / D1;
    float z2 = y2 / D2;
    float z3 = y3 / D3;
    float z4 = y4 / D4;
    float z5 = y5 / D5;
    float z6 = y6 / D6;

    // Step 4: Backward substitution: Solve L^T x = z
    xAng[2] = z6;
    xAng[1] = z5 - L65 * xAng[2];
    xAng[0] = z4 - L54 * xAng[1] - L64 * xAng[2];
    xLin[2] = z3 - L43 * xAng[0] - L53 * xAng[1] - L63 * xAng[2];
    xLin[1] = z2 - L32 * xLin[2] - L42 * xAng[0] - L52 * xAng[1] - L62 * xAng[2];
    xLin[0] = z1 - L21 * xLin[1] - L31 * xLin[2] - L41 * xAng[0] - L51 * xAng[1] - L61 * xAng[2];
}

} // namespace avbd

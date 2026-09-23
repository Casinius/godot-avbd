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
 * Port note: an earlier revision of this file contained `using namespace std;`,
 * which leaked every standard-library name into any translation unit that
 * included it. The handful of places that relied on it now name `std::`
 * explicitly.
 */

#pragma once

#include "Eigen/Core"
#include "Eigen/Geometry"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <cmath>
#include <cstddef>

namespace avbd {

// -----------------------------------------------------------------------------
// Math types
//
// Plain aggregates: no constructors, no invariants to break, so they stay
// trivially copyable and cheap to pass by value. Component access goes through
// `operator[]`, which spells the field out rather than aliasing the struct's
// address as an array - that trick is undefined behaviour, and it was the only
// reason this header needed the C-style casts it used to contain. With a
// constant index the comparison folds away entirely.
// -----------------------------------------------------------------------------

// struct float2
// {
//     float x, y;

//     float &operator[](std::size_t i) noexcept { return i == 0 ? x : y; }
//     const float &operator[](std::size_t i) const noexcept { return i == 0 ? x
//     : y; }
// };

// struct float3
// {
//     float x, y, z;

//     float &operator[](std::size_t i) noexcept { return i == 0 ? x : (i == 1 ?
//     y : z); } const float &operator[](std::size_t i) const noexcept { return
//     i == 0 ? x : (i == 1 ? y : z); }
// };

// struct quat
// {
//     float x, y, z, w;

//     float &operator[](std::size_t i) noexcept { return i == 0 ? x : (i == 1 ?
//     y : (i == 2 ? z : w)); } const float &operator[](std::size_t i) const
//     noexcept { return i == 0 ? x : (i == 1 ? y : (i == 2 ? z : w)); }
// };

// struct float2x2
// {
//     float2 row[2];

//     float2 &operator[](std::size_t i) noexcept { return row[i]; }
//     const float2 &operator[](std::size_t i) const noexcept { return row[i]; }

//     [[nodiscard]] float2 col(std::size_t i) const noexcept { return
//     float2{row[0,i], row[1,i]}; }
// };

// struct float3x3
// {
//     float3 row[3];

//     float3 &operator[](std::size_t i) noexcept { return row[i]; }
//     const float3 &operator[](std::size_t i) const noexcept { return row[i]; }

//     [[nodiscard]] float3 col(std::size_t i) const noexcept { return
//     float3{row[0,i], row[1,i], row[2,i]}; }
// };

using float3x3 = Eigen::Matrix<float, 3, 3>;
using float3 = Eigen::Vector3<float>;

using quat = Eigen::Quaternion<float>;
[[nodiscard]] inline float3x3 diagonal(float m00, float m11, float m22) noexcept
{
    float3x3 m;
    m << m00, 0.0f, 0.0f,
         0.0f, m11, 0.0f,
         0.0f, 0.0f, m22;
    return m;
}
// //
// -----------------------------------------------------------------------------
// // quat operators
// //
// -----------------------------------------------------------------------------

// [[nodiscard]] inline quat operator*(quat a, float b) noexcept
// {
//     return {a.x() * b, a.y() * b, a.z() * b, a.w() * b};
// }

// [[nodiscard]] inline quat operator/(quat a, float b) noexcept
// {
//     return {a.x() / b, a.y() / b, a.z() / b, a.w() / b};
// }

// [[nodiscard]] inline quat operator*(quat a, quat b) noexcept
// {
//     return {
//         a.w() * b.x() + a.x() * b.w() + a.y() * b.z() - a.z() * b.y(),
//         a.w() * b.y() - a.x() * b.z() + a.y() * b.w() + a.z() * b.x(),
//         a.w() * b.z() + a.x() * b.y() - a.y() * b.x() + a.z() * b.w(),
//         a.w() * b.w() - a.x() * b.x() - a.y() * b.y() - a.z() * b.z()};
// }

// [[nodiscard]] inline quat operator+(quat a, quat b) noexcept
// {
//     return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z(), a.w() + b.w()};
// }
[[nodiscard]] inline float3 operator-(quat a, quat b) noexcept {
  const quat d = a * b.inverse();
  return float3{d.x(), d.y(), d.z()} * 2.0f;
}

// Solve the symmetric 6x6 system [aLin  aCross^T; aCross  aAng] x = [bLin;
// bAng] using Eigen::LDLT. The matrix is stored as its lower triangle; we
// reconstruct a full 6x6 matrix for Eigen.
inline void solve(float3x3 aLin, float3x3 aAng, float3x3 aCross, float3 bLin,
                  float3 bAng, float3 &xLin, float3 &xAng) noexcept {

  // Build 6x6 symmetric matrix [aLin  aCross^T; aCross  aAng] from custom
  // structs
  Eigen::MatrixXf A(6, 6);
  // Top-left: aLin (lower triangle only)
  A.setZero();
  // 把 aLin / aAng 填到对应块，aCross 填到对应块，并且对称化
  A.block<3, 3>(0, 0) = aLin;
  A.block<3, 3>(3, 3) = aAng;
  A.block<3, 3>(0, 3) = aCross.transpose();
  A.block<3, 3>(3, 0) = aCross;
  // Stack bLin and bAng into 6-vector
  Eigen::MatrixXf b(6, 1);
  b << bLin[0], bLin[1], bLin[2], bAng[0], bAng[1], bAng[2];

  // Solve using LDLT decomposition
  Eigen::LDLT<Eigen::MatrixXf> solver(A);
  Eigen::MatrixXf x = solver.solve(b);

  // Extract solution into output custom structs (use parentheses for Eigen)
  xLin[0] = x(0, 0);
  xLin[1] = x(1, 0);
  xLin[2] = x(2, 0);
  xAng[0] = x(3, 0);
  xAng[1] = x(4, 0);
  xAng[2] = x(5, 0);
}

} // namespace avbd

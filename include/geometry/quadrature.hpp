/// @file quadrature.hpp
///
/// the class related to quadrature.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_GEOMETRY_QUADRATURE_HPP_
#define INCLUDE_GEOMETRY_QUADRATURE_HPP_

#include "defs.hpp"
#include <array>

using namespace std;

template <int kPoint>
struct LineQuad;

template <>
struct LineQuad<1> {
  array<Real, 1> x{0.0};
  array<Real, 1> w{2.0};
};
template <>
struct LineQuad<2> {
  array<Real, 2> x{-0.577350269189626, 0.577350269189626};
  array<Real, 2> w{1.0, 1.0};
};
template <>
struct LineQuad<3> {
  array<Real, 3> x{-0.774596669241483, +0.774596669241483, 0.0};
  array<Real, 3> w{+0.555555555555556, +0.555555555555556, 0.888888888888889};
};
template <>
struct LineQuad<4> {
  array<Real, 4> x{-0.8611363115940520, +0.8611363115940520, 
                   -0.3399810435848560, +0.3399810435848560};
  array<Real, 4> w{+0.3478548451374530, +0.3478548451374530,
                   +0.6521451548625460, +0.6521451548625460};
};

template <int kOrder>
struct TriQuad;

template <>
struct TriQuad<2> {
  array<Real, 3> a{0.0, 0.5, 0.5};
  array<Real, 3> b{0.5, 0.0, 0.5};
  array<Real, 3> c{0.5, 0.5, 0.0};
  array<Real, 3> w{0.3333333333333333, 0.3333333333333333, 0.3333333333333333};
};

template <>
struct TriQuad<3> {
  array<Real, 4> a{0.3333333333333333, 0.6, 0.2, 0.2};
  array<Real, 4> b{0.3333333333333333, 0.2, 0.2, 0.6};
  array<Real, 4> c{0.3333333333333333, 0.2, 0.6, 0.2};
  array<Real, 4> w{-0.5624999999999998, 0.5208333333333332,
                                        0.5208333333333332,
                                        0.5208333333333332};
};

template <int kOrder>
struct QuaQuad;

template <>
struct QuaQuad<2> {
  array<Real, 4> a{+0.577350269189626, -0.577350269189626, -0.577350269189626, +0.577350269189626};
  array<Real, 4> b{+0.577350269189626, +0.577350269189626, -0.577350269189626, -0.577350269189626};
  array<Real, 4> w{1.0, 1.0, 1.0, 1.0};
};

template <>
struct QuaQuad<3> {
  array<Real, 4> a{+0.577350269189626, -0.577350269189626, -0.577350269189626, +0.577350269189626};
  array<Real, 4> b{+0.577350269189626, +0.577350269189626, -0.577350269189626, -0.577350269189626};
  array<Real, 4> w{1.0, 1.0, 1.0, 1.0};
};

#endif  //  INCLUDE_GEOMETRY_QUADRATURE_HPP_

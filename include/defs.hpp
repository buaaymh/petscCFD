/// @file defs.h
///
/// Global definitions, structures and macros.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_DSFS_HPP_
#define INCLUDE_DSFS_HPP_

#include <cmath>
#include <cfloat>
#include <Eigen/Dense>
#include <petscdmplex.h>

// #include "solver.hpp"

#define DIM 2   /* Geometric dimension */

  /// Kind of equations solved
  enum class Equations { Linear, Euler, NavierStokes };

  /// Kind of time-stepping
  enum class TimeStepping { RK3 };

  /// Kind of edge and boundary
  enum class BdCondType { Interior, Periodic, InFlow, OutFlow,
                          FarField, InviscWall, Symmetry };

// floating point type (SGLPREC=single precision, otherwise double) ***********

// #define SGLPREC

#ifdef SGLPREC
  #define PETSC_USE_REAL_SINGLE
  typedef float  Real;        /**< floating-point value */
  #define Abs    fabsf
  #define Sqrt   sqrtf
  #define Sin    sinf
  #define Cos    cosf
  #define Pow    powf
  #define Tan    tanf
#else
  #define PETSC_USE_REAL_DOUBLE
  typedef double Real;        /**< floating-point value */
  #define Abs    std::fabs
  #define Sqrt   std::sqrt
  #define Sin    std::sin
  #define Cos    std::cos
  #define Pow    std::pow
  #define Tan    std::tan
#endif

// general constants **********************************************************

#define PI                         3.1415926535897932
#define Gamma                      1.4
#define GammaPlusOne               Gamma+1
#define GammaMinusOne              Gamma-1
#define OneOverGamma               1/Gamma
#define OneOverGammaMinusOne       1/GammaMinusOne

// declare template class macros ***********************************************
template<int kOrder, class Physics>
struct BndConds;
template<int kOrder, class Physics>
struct EdgeGroup;
template<int kOrder, class Physics>
class Solver;
template <int kOrder>
class Edge;
template <int kOrder>
class Cell;
template <int kOrder>
class Mesh;

// static function class macros ***********************************************

template <int kOrder>
struct Dp {
  static constexpr Real Factorial(int p) {
    int fac = 1;
    for (int i = 1; i <= p; ++i) { fac *= i; }
    return Real(fac);
  }
  static void GetDpArrayInterior(Real distance, Real* dp) {
    for (int i = 0; i <= kOrder; ++i) {
      dp[i] = Pow(distance, 2*i-1) / Pow(Factorial(i), 2);
    }
  }
  static void GetDpArrayP0(Real distance, Real* dp) {
    dp[0] = 1 / distance;
    for (int i = 1; i <= kOrder; ++i) { dp[i] = 0; }
  }
  static void GetDpArraySymmetry(Real distance, Real* dp) {
    for (int i = 0; i <= kOrder; ++i) {
      dp[i] = Pow(distance, 2*i-1) * Pow(-1,i) / Pow(Factorial(i), 2);
    }
  }
};

template <int nEqual, int kOrder>
struct Algebra {
  using State = Eigen::Matrix<Real, nEqual, 1>;
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr State Grad(const Real* basis, const Real* valCoefs) {
    State state = State::Zero();
    int iter{0};
    for (int i = 0; i < nEqual; ++i) {
      for (int j = 0; j < nCoef; ++j) {
        state(i) += valCoefs[iter++] * basis[j];
      }
    }
    return state;
  }
  static constexpr void Add(int n, Real* a, const Real* b) {
    for (int i = 0; i < n; ++i) { a[i] += b[i]; }
  }
  static constexpr void Sub(int n, Real* a, const Real* b) {
    for (int i = 0; i < n; ++i) { a[i] -= b[i]; }
  }
};

#endif // INCLUDE_DSFS_HPP_
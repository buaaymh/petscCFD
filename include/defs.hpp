/// @file defs.hpp
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

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <Eigen/Dense>
#include <petscdmplex.h>
#include <petscts.h>

#define DIM 2   /* Geometric dimension */

  /// Kind of equations solved
  enum class Equations { Linear, Euler, NavierStokes };

  /// Kind of time-stepping
  enum class TimeStepping { RK3 };

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

// function macros ************************************************************

static constexpr Real Sign(Real a)
{
  if (a < 0) { return -1; }
  else if (a > 0) { return 1;}
  else { return 0; }
}

static constexpr Real Minmod(Real a, Real b)
{
  if (a * b <= 0) { return 0; }
  else { return Sign(a) * std::min(Abs(a), Abs(b)); }
}

#endif // INCLUDE_DSFS_HPP_
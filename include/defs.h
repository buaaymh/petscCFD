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

#ifndef INCLUDE_DSFS_H_
#define INCLUDE_DSFS_H_
// STL
#include <array>
#include <cmath>
#include <cfloat>
#include <functional>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
// Eigen3
#include <Eigen/Dense>
// PETSc 
#include <petscdmplex.h>
#include <petscsf.h>

namespace cfd {

using std::array;
using std::cout;
using std::endl;
using std::function;
using std::make_unique;
using std::max;
using std::min;
using std::set;
using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

using Eigen::Matrix;
using Eigen::Dynamic;

#define DIM 2   /* Geometric dimension */

  /// Kind of equations solved
  enum class Equations { Linear, Euler, NavierStokes };

  /// Kind of time-stepping
  enum class TimeStepping { Global, Local };

  /// Kind of edge and boundary
  enum class BdCondType { Interior, Periodic, InFlow, OutFlow,
                          FarField, InviscWall, Symmetry };

// floating point type (SGLPREC=single precision, otherwise double) ***********

#ifdef SGLPREC
  #define PETSC_USE_REAL_SINGLE
  typedef float  Real;        /**< floating-point value */
  #define EPS    FLT_EPSILON  /**< tolerance for a small number */
  #define INF    FLT_MAX      /**< tolerance for a large number */
  #define Abs    fabsf
  #define Sqrt   sqrtf
  #define Sin    sinf
  #define Cos    cosf
  #define Pow    powf
  #define Tan    tanf
#else
  #define PETSC_USE_REAL_DOUBLE
  typedef double Real;        /**< floating-point value */
  #define EPS    DBL_EPSILON  /**< tolerance for a small number */
  #define INF    DBL_MAX      /**< tolerance for a large number */
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
#define GamP1                      2.4
#define GamM1                      0.4
#define One_Gamma                  0.7142857142857143
#define One_GamM1                  2.5
#define Gam_GamM1                  3.5
#define EPS_Entr                   0.05

// declare template class macros ***********************************************
template<int kOrder, class Physics>
struct Group;
template<int kOrder, class Physics>
class Solver;
template <int kOrder, class Physics>
class SpaceDiscr;
template <int kOrder>
class Edge;
template <int kOrder>
class Cell;
template <int kOrder>
class Mesh;

// declare some simple class **************************************************

typedef Eigen::Matrix<Real, 2, 1> Node;

typedef Eigen::Matrix<Real, Dynamic, Dynamic> Array;

template <int kOrder>
struct cmp {
  bool operator()(Edge<kOrder>* a, Edge<kOrder>* b) const {
    return a->I() < b->I();
  }
};
template <int kOrder>
using EdgeSet = set<Edge<kOrder>*, cmp<kOrder>>;

template<int kOrder, class Physics>
using EdgeGroups = unordered_map<int, unique_ptr<Group<kOrder,Physics>>>;

class BndConds {
 public:
  // For periodic boundary ****************************************************
  Real lower[2];
  Real upper[2];
  // For far-field boundary ****************************************************
  Real* refVal;
  // For inflow boundary ******************************************************
  std::function<void(Real, const Real*, Real*)> inflow;

  BndConds();
  ~BndConds();
};

// static function class macros ***********************************************
static constexpr Real Factorial(int p) {
  int fac = 1;
  for (int i = 1; i <= p; ++i) { fac *= i; }
  return static_cast<Real>(fac);
}
static constexpr void InteriorDp(int k, Real distance, Real* dp) {
  for (int i = 0; i <= k; ++i) {
    dp[i] = Pow(distance, 2*i-1) / Pow(Factorial(i), 2);
  }
}
static constexpr void WithoutDerivative(int k, Real distance, Real* dp) {
  dp[0] = 1 / distance;
  for (int i = 1; i <= k; ++i) { dp[i] = 0; }
}
static constexpr void Symmetry(int k, Real distance, Real* dp) {
  for (int i = 0; i <= k; ++i) {
    dp[i] = Pow(distance, 2*i-1) * Pow(-1,i) / Pow(Factorial(i), 2);
  }
}
static constexpr Real EntropyCorr(Real z, Real d)
{
  if (z > d) return z;
  else return 0.5*(z*z+d*d)/d;
}

static constexpr Real Venkat(Real d1, Real d2, Real eps2)
{
  Real num = (d1 * d1 + eps2) * d2 + 2.0 * d2 * d2 * d1;
  Real den = d2 * (d1 * d1 + 2.0 * d2 * d2 + d1 * d2 + eps2);
  return num / den;
}

}  // cfd

#endif // INCLUDE_DSFS_H_
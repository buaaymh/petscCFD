/// @file physModel.hpp
///
/// Definition of the class related to physical models.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_PHYSMODEL_HPP_
#define INCLUDE_PHYSMODEL_HPP_

#include "defs.hpp"

using namespace std;

struct Linear {
  Linear() = default;
  static constexpr int nEqual = 1;
  using Flux = Eigen::Matrix<Real, nEqual, 1>;
  using State = Flux;
  static unordered_map<string, int> CreateFieldDiscription() {
    unordered_map<string, int>  field_desc;
    field_desc.emplace("U", 1);
    return field_desc;
  }
  static Flux Riemann(int dim, const Real* coord, const Real* normal,
                                const Real* U_l, const Real* U_r) {
    Flux F;
    if (normal[0] > 0) { F(0) = U_l[0] * normal[0]; }
    else { F(0) = U_r[0] * normal[0]; }
    return F;
  }
};
struct Euler
{
  Euler() = default;
  static constexpr int nEqual = 4;
  using Flux = Eigen::Matrix<Real, nEqual, 1>;
  using State = Flux;
  static unordered_map<string, int> CreateFieldDiscription() {
    unordered_map<string, int>  field_desc;
    field_desc.emplace("Density", 1);
    field_desc.emplace("U", 1);
    field_desc.emplace("V", 1);
    field_desc.emplace("Pressure", 1);
    return field_desc;
  }
};

#endif // INCLUDE_PHYSMODEL_HPP_

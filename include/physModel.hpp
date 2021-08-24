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

struct Linear
{
  Linear() = default;
  static constexpr int nEqual = 1;
  static unordered_map<string, int> CreateFieldDiscription() {
    unordered_map<string, int>  field_desc;
    field_desc.emplace("Density", 1);
    return field_desc;
  }
};
struct Euler
{
  Euler() = default;
  static constexpr int nEqual = 4;
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

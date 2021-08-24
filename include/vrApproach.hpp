/// @file vrApproach.hpp
///
/// Definition of the class related to variational reconstruction approach.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_VRAPPROACH_HPP_
#define INCLUDE_VRAPPROACH_HPP_

#include "defs.hpp"
#include "geometry/mesh.hpp"

using namespace std;

template <int kOrder, class Physics>
class VrApproach
{
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
 public:
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  using Column = Eigen::Matrix<Real, nCoef, 1>;
  using MeshType = Mesh<2, kOrder>;
  using BndCondsType = BndConds<MeshType, Physics>;
  typedef struct T_VrBlock {
    Matrix C_mat;
    Column b_sub;
  } VrBlock;

  vector<Matrix> A_inv;
  vector<Matrix> B_mat;
  vector<Column> b_col;
  vector<VrBlock> block;

  // functions
  VrApproach() = default;
  void BuildBlock(const MeshType& mesh, const BndCondsType& bndcond);

 private:
};

#endif // INCLUDE_VRAPPROACH_HPP_

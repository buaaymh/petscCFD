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
  using MeshType = Mesh<kOrder>;
  using BndCondsType = BndConds<MeshType, Physics>;
  typedef struct T_VrBlock {
    Matrix C_mat;
    Column b_sub;
  } VrBlock;
  DM            dmCoef;
  PetscSF       sfCoef;
  void SetCoefLayout(DM dm) {
    PetscSection    section;
    int             cStart, cEnd;
    int             nroots, nleaves;
    
    DMClone(dm, &dmCoef);
    PetscSectionCreate(PetscObjectComm((PetscObject)dm), &section);
    DMPlexGetHeightStratum(dmCoef, 0, &cStart, &cEnd);
    PetscSectionSetChart(section, cStart, cEnd);
    for (int c = cStart; c < cEnd; ++c)
      PetscSectionSetDof(section, c, nCoef * Physics::nEqual);
    PetscSectionSetUp(section);
    DMSetLocalSection(dmCoef, section);
    PetscSectionDestroy(&section);
    DMGetSectionSF(dmCoef, &sfCoef);
    /* Build the ghosted start forest for data */
    DMGetSectionSF(dmCoef, &sfCoef);
    PetscSFGetGraph(sfCoef, &nroots, &nleaves, nullptr, nullptr);
    int selected[nleaves-nroots];
    for (int i = nroots; i < nleaves; ++i) { selected[i-nroots] = i; }
    PetscSFCreateEmbeddedLeafSF(sfCoef, nleaves-nroots, selected, &sfCoef);
  }

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

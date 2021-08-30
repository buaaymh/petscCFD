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

#include "defs.h"

namespace cfd {

template <int kOrder, class Physics>
class VrApproach {
 public:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */

  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  using Column = Eigen::Matrix<Real, nCoef, 1>;
  using EqualCol = Eigen::Matrix<Real, nCoef, Physics::nEqual>;
  using MeshType = Mesh<kOrder>;
  using Cell = typename MeshType::CellType;
  using FuncTable = typename Cell::BasisF;
  using Manager = EdgeManager<kOrder, Physics>;
  using Set = typename EdgeGroup<kOrder, Physics>::EdgeSet;
  using Coefs = Eigen::Matrix<Real, nCoef, Dynamic>;

  DM            dmCoef;
  PetscSF       sfCoef;
  Coefs         coefs;

  struct VrBlock {
    VrBlock() : C_mat(Matrix::Zero()), b_sub(Column::Zero()) {}
    Matrix C_mat;
    Column b_sub;
  };
  vector<Matrix> A_inv;
  vector<Matrix> B_mat;
  vector<EqualCol> b_col;
  vector<VrBlock> block;
  VrApproach() = default;
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
    // Initialize coefficients
    coefs.resize(nCoef, cEnd * Physics::nEqual);
    coefs.setZero();
  }
  void AllocatorMats(const MeshType& mesh) {
    int n_cell = mesh.NumLocalCells(), n_edge = mesh.CountEdges();
    A_inv = vector<Matrix>(n_cell, Matrix::Zero());
    B_mat = vector<Matrix>(n_edge, Matrix::Zero());
    b_col = vector<EqualCol>(n_cell, EqualCol::Zero());
    block = vector<VrBlock>(mesh.neighbor.size(), VrBlock());
  }
  void CalculateAinvs(const MeshType& mesh, const Manager& edgeManager) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      for (int j = mesh.offset[i]; j < mesh.offset[i+1]; ++j) {
        auto e = mesh.edge[mesh.interface[j]].get();
        Matrix temp = Matrix::Zero();
        Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1];
        if (mesh.neighbor[j] >= 0) { // Interiod, periodic and outflow boundary
          InteriorDp(kOrder, distance, dp);
        } else {
          BdCondType type = static_cast<BdCondType>(-mesh.neighbor[j]);
          if (type == BdCondType::InFlow ||
              type == BdCondType::FarField ||
              type == BdCondType::InviscWall) {
            WithoutDerivative(kOrder, distance, dp);
          } else if (type == BdCondType::Symmetry) {
            Symmetry(kOrder, distance, dp);
          }
        }
        e->Integrate([&](const Node& node) {
          return GetMatAt(node.data(), *c, *c, normal, dp);
        }, &temp);
        A_inv[i] += temp;
        e->Integrate([&](const Node& node) {
          return GetVecAt(*c, node.data(), distance); }, &(block[j].b_sub));
      }
      A_inv[i] = A_inv[i].inverse();
    }
  }
  void CalculateBlockC(const MeshType& mesh) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      for (int j = mesh.offset[i]; j < mesh.offset[i+1]; ++j) {
        if (mesh.neighbor[j] < i) {
          block[j].C_mat = A_inv[i] * B_mat[mesh.interface[j]];
        } else {
          block[j].C_mat = A_inv[i] * B_mat[mesh.interface[j]].transpose();
        }
      }
    }
  }
  static Matrix GetMatAt(const Real* coord, const Cell& a, const Cell& b,
                         const Real* normal, const Real* dp) {
    FuncTable i = a.GetFuncTable(coord, normal);
    FuncTable j = b.GetFuncTable(coord, normal);
    Matrix mat = Matrix::Zero();
    for (int m = 0; m != nCoef; ++m) {
      for (int n = 0; n != nCoef; ++n) {
        for (int k = 0; k != kOrder+1; ++k) mat(m,n) += dp[k] * i(n,k) * j(m,k);
      }
    }
    if (a.I() > b.I()) mat.transposeInPlace();
    return mat;
  }

 private:
  static Column GetVecAt(const Cell& a, const Real* coord, Real distance) {
    return a.Functions(coord) / distance;
  }
};

}  // cfd

#endif // INCLUDE_VRAPPROACH_HPP_

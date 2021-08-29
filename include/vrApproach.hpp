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

#include "bndConds.hpp"
#include "defs.hpp"
#include "geometry/mesh.hpp"
#include <iostream>

namespace cfd {

using std::vector;

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

  DM            dmCoef;
  PetscSF       sfCoef;

  struct VrBlock {
    VrBlock() : C_mat(Matrix::Zero()), b_sub(Column::Zero()) {}
    Matrix C_mat;
    Column b_sub;
  };

  vector<int> offset;
  vector<Matrix> A_inv;
  vector<Matrix> B_mat;
  vector<EqualCol> b_col;
  vector<VrBlock> block;

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
  void AllocatorMats(const MeshType& mesh) {
    int n_cell = mesh.NumLocalCells(), n_edge = mesh.CountEdges();
    A_inv = vector<Matrix>(n_cell, Matrix::Zero());
    B_mat = vector<Matrix>(n_edge, Matrix::Zero());
    b_col = vector<EqualCol>(n_cell, EqualCol::Zero());
    offset.reserve(n_cell);
    offset.emplace_back(0);
    for (int i = 1; i < n_cell; ++i) {
      offset.emplace_back(offset[i-1] + mesh.cell[i-1]->nCorner());
    }
    int n = offset[n_cell-1] + mesh.cell[n_cell-1]->nCorner();
    block = vector<VrBlock>(n, VrBlock());
  }
  void CalculateAinvs(const MeshType& mesh, const Manager& edgeManager) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      for (int j = 0; j < c->nCorner(); ++j) {
        auto e = mesh.edge[c->Edge(j)].get();
        Matrix temp = Matrix::Zero();
        Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1];
        if (c->Adjc(j) >= 0) { // Interiod, Periodic boundary
          Dp<kOrder>::Interior(distance, dp);
        } else {
          edgeManager.bdGroup.at(-c->Adjc(j))->GetDpArray(distance, dp);
        }
        e->Integrate([&](const Node& node) {
          return GetMatAt(node.data(), *c, *c, normal, dp);
        }, &temp);
        A_inv[i] += temp;
        e->Integrate([&](const Node& node) {
          return GetVecAt(*c, node.data(), distance);
        }, &(block[offset[i]+j].b_sub));
      }
      A_inv[i] = A_inv[i].inverse();
    }
  }
  void CalculateBlockC(const MeshType& mesh) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      for (int j = 0; j < c->nCorner(); ++j) {
        if (c->Adjc(j) < i) {
          block[offset[i]+j].C_mat = A_inv[i] * B_mat[c->Edge(j)];
        } else {
          block[offset[i]+j].C_mat = A_inv[i] * B_mat[c->Edge(j)].transpose();
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

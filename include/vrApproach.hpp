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
  static constexpr int nEqual = Physics::nEqual;
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  using Column = Eigen::Matrix<Real, nCoef, 1>;
  using EqualCol = Eigen::Matrix<Real, nCoef, nEqual>;
  using MeshType = Mesh<kOrder>;
  using Cell = typename MeshType::CellType;
  using FuncTable = typename Cell::BasisF;
  using Manager = EdgeManager<kOrder, Physics>;
  using Set = typename EdgeGroup<kOrder, Physics>::EdgeSet;
  using ConVar = typename Eigen::Matrix<Real, nEqual, Dynamic>;

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
    PetscSectionCreate(PetscObjectComm((PetscObject)dmCoef), &section);
    DMPlexGetHeightStratum(dmCoef, 0, &cStart, &cEnd);
    PetscSectionSetChart(section, cStart, cEnd);
    for (int c = cStart; c < cEnd; ++c)
      PetscSectionSetDof(section, c, nCoef * nEqual);
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
    coefs.resize(nCoef, cEnd * nEqual);
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
      Matrix a_mat = Matrix::Zero();
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
        a_mat += temp;
        e->Integrate([&](const Node& node) {
          return GetVecAt(*c, node.data(), distance); }, &(block[j].b_sub));
      }
      A_inv[i] = a_mat.inverse();
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
  void UpdateBcols(const MeshType& mesh, const ConVar& cv, const Manager& edgeManager) {
    const auto& offset = mesh.offset;
    const auto& interface = mesh.interface;
    const auto& neighbor = mesh.neighbor;
    for (int i = 0; i < b_col.size(); ++i) {
      b_col[i].setZero();
      for (int j = offset[i]; j < offset[i+1]; ++j) {
        if (neighbor[j] >= 0) {
          auto cv_d = cv.col(neighbor[j]) - cv.col(i);
          b_col[i] += block[j].b_sub * cv_d.transpose();
        } else {
        }
      }
      b_col[i] = A_inv[i] * b_col[i];
    }
  }
  void UpdateCoefs(const MeshType& mesh, const Manager& edgeManager) {
    const auto& offset = mesh.offset;
    const auto& interface = mesh.interface;
    const auto& neighbor = mesh.neighbor;
    for (int k = 0; k < 8; ++k) {
      for (int i = 0; i < b_col.size(); ++i) {
        coefs.block<nCoef, nEqual>(0, i*nEqual) *= -0.3;
        EqualCol temp = EqualCol::Zero();
        for (int j = offset[i]; j < offset[i+1]; ++j) {
          if (neighbor[j] >= 0) {
            temp += block[j].C_mat * coefs.block<nCoef, nEqual>(0, neighbor[j]*nEqual);
          } else {
          }
        }
        temp += b_col[i];
        coefs.block<nCoef, nEqual>(0, i*nEqual) += temp * 1.3;
      }
      PetscSFScatterBegin(sfCoef, MPIU_REAL, coefs.data(), coefs.data());
      PetscSFScatterEnd(sfCoef, MPIU_REAL, coefs.data(), coefs.data());
    }
  }
 private:
  static Column GetVecAt(const Cell& a, const Real* coord, Real distance) {
    return a.Functions(coord) / distance;
  }
};

}  // cfd

#endif // INCLUDE_VRAPPROACH_HPP_

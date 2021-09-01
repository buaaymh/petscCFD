/// @file spaceDiscr.hpp
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

#ifndef INCLUDE_SPACEDISCR_HPP_
#define INCLUDE_SPACEDISCR_HPP_

#include "defs.h"

namespace cfd {

template <int kOrder, class Physics>
class SpaceDiscr
{
 public:
  // Constants:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr int nEqual = Physics::nEqual;
  // Types:
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  using Column = Eigen::Matrix<Real, nCoef, 1>;
  using EqualCol = Eigen::Matrix<Real, nCoef, nEqual>;
  using EdgeGroups = typename Solver<kOrder, Physics>::EdgeGroups;

  DM            dmCoef;
  PetscSF       sfCoef;
  Array         conVar;
  Array         coefs;

  struct VrBlock {
    VrBlock() : C_mat(Matrix::Zero()), b_sub(Column::Zero()) {}
    Matrix C_mat;
    Column b_sub;
  };
  vector<Matrix> A_inv;
  vector<Matrix> B_mat;
  vector<EqualCol> b_col;
  vector<VrBlock> block;
  // Functions:
  SpaceDiscr() = default;
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
  void AllocatorMats(const Mesh<kOrder>& mesh) {
    int n_cell = mesh.NumLocalCells(), n_edge = mesh.CountEdges();
    A_inv = vector<Matrix>(n_cell, Matrix::Zero());
    B_mat = vector<Matrix>(n_edge, Matrix::Zero());
    b_col = vector<EqualCol>(n_cell, EqualCol::Zero());
    block = vector<VrBlock>(mesh.adjc_csr.size(), VrBlock());
  }
  void CalculateBmats(const EdgeGroups& faceGroups) {
    for (auto& [type, group] : faceGroups) {
      if (BdCondType(type) == BdCondType::Interior ||
          BdCondType(type) == BdCondType::Periodic) {
        for (auto& e : group->edge) {
          B_mat[e->I()] = group->CalculateMat(e, e->left, e->right);
        }
      } else {
        for (auto& e : group->edge) {
          B_mat[e->I()] = group->CalculateMat(e, e->left, e->left);
        }
      }
    }
  }
  void CalculateAinvs(const Mesh<kOrder>& mesh, const EdgeGroups& faceGroups) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      Matrix a_mat = Matrix::Zero();
      for (int j = mesh.offset[i]; j < mesh.offset[i+1]; ++j) {
        auto e = mesh.edge[mesh.edge_csr[j]].get();
        Matrix temp = Matrix::Zero();
        Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1];
        if (mesh.adjc_csr[j] >= 0) { // Interiod, periodic and outflow boundary
          a_mat += faceGroups.at(0)->CalculateMat(e, c, c);
        } else {
          a_mat += faceGroups.at(-mesh.adjc_csr[j])->CalculateMat(e, c, c);
        }
        e->Integrate([&](const Node& node) {
          return GetVecAt(*c, node.data(), distance); }, &(block[j].b_sub));
      }
      A_inv[i] = a_mat.inverse();
    }
  }
  void CalculateBlockC(const Mesh<kOrder>& mesh) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      for (int j = mesh.offset[i]; j < mesh.offset[i+1]; ++j) {
        if (mesh.adjc_csr[j] < i) {
          block[j].C_mat = A_inv[i] * B_mat[mesh.edge_csr[j]];
        } else {
          block[j].C_mat = A_inv[i] * B_mat[mesh.edge_csr[j]].transpose();
        }
      }
    }
  }
  static constexpr Matrix GetMatAt(const Real* coord, const Cell<kOrder>& a, const Cell<kOrder>& b,
                         const Real* normal, const Real* dp) {
    auto i = a.GetFuncTable(coord, normal);
    auto j = b.GetFuncTable(coord, normal);
    Matrix mat = Matrix::Zero();
    for (int m = 0; m != nCoef; ++m) {
      for (int n = 0; n != nCoef; ++n) {
        for (int k = 0; k != kOrder+1; ++k) mat(m,n) += dp[k] * i(n,k) * j(m,k);
      }
    }
    if (a.I() > b.I()) mat.transposeInPlace();
    return mat;
  }
  void UpdateBcols(const Mesh<kOrder>& mesh, const Array& conVar) {
    const auto& offset = mesh.offset;
    const auto& adjc_csr = mesh.adjc_csr;
    for (int i = 0; i < b_col.size(); ++i) {
      b_col[i].setZero();
      for (int j = offset[i]; j < offset[i+1]; ++j) {
        if (adjc_csr[j] >= 0) {
          auto conVar_d = conVar.col(adjc_csr[j]) - conVar.col(i);
          b_col[i] += block[j].b_sub * conVar_d.transpose();
        } else {
        }
      }
      b_col[i] = A_inv[i] * b_col[i];
    }
  }
  void UpdateCoefs(const Mesh<kOrder>& mesh) {
    const auto& offset = mesh.offset;
    const auto& adjc_csr = mesh.adjc_csr;
    for (int k = 0; k < 8; ++k) {
      for (int i = 0; i < b_col.size(); ++i) {
        coefs.block<nCoef, nEqual>(0, i*nEqual) *= -0.3;
        EqualCol temp = EqualCol::Zero();
        for (int j = offset[i]; j < offset[i+1]; ++j) {
          if (adjc_csr[j] >= 0) {
            temp += block[j].C_mat * coefs.block<nCoef, nEqual>(0, adjc_csr[j]*nEqual);
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
  static Column GetVecAt(const Cell<kOrder>& a, const Real* coord, Real distance) {
    return a.Functions(coord) / distance;
  }
};

}  // cfd

#endif // INCLUDE_SPACEDISCR_HPP_

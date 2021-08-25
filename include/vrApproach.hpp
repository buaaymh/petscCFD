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
#include "bndConds.hpp"

using namespace std;

template <int kOrder, class Physics>
class VrApproach
{
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
 public:
  using Matrix = Eigen::Matrix<float, nCoef, nCoef>;
  using Column = Eigen::Matrix<float, nCoef, 1>;
  using MeshType = Mesh<kOrder>;
  using Cell = typename MeshType::CellType;
  using FuncTable = typename Cell::BasisF;
  using Set = typename MeshType::EdgeSet;
  using BndCondsType = BndConds<MeshType, Physics>;
  struct VrBlock {
    VrBlock() : C_mat(Matrix::Zero()), b_sub(Column::Zero()) {}
    Matrix C_mat;
    Column b_sub;
  };
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
  void AllocatorMats(const MeshType& mesh) {
    int n_cell = mesh.NumLocalCells(), n_edge = mesh.CountEdges();
    A_inv = vector<Matrix>(n_cell, Matrix::Zero());
    B_mat = vector<Matrix>(n_edge, Matrix::Zero());
    b_col = vector<Column>(n_cell, Column::Zero());
    offset.reserve(n_cell);
    offset.emplace_back(0);
    for (int i = 1; i < n_cell; ++i) {
      offset.emplace_back(offset[i-1] + mesh.cell[i-1]->nCorner());
    }
    int n = offset[n_cell-1] + mesh.cell[n_cell-1]->nCorner();
    block = vector<VrBlock>(n, VrBlock());
  }
  void CalculateBmats(const Set& edges, BdCondType type) {
    if (type == BdCondType::Interior ||
        type == BdCondType::Periodic ||
        type == BdCondType::OutFlow) {
      for (auto& e : edges) {
        Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1];
        for (int i = 0; i <= kOrder; ++i){
          dp[i] = Pow(distance, 2*i-1)/ Pow(Factorial(i), 2);
        }
        e->Integrate([&](const Node& node) {
          return GetMatAt(node.data(), *(e->left), *(e->right), normal, dp);
        }, &B_mat[e->I()]);
        // cout << B_mat[e->I()] << endl << endl;
      }
    }
  }
  void CalculateAinvs(const MeshType& mesh) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      for (int j = 0; j < c->nCorner(); ++j) {
        auto e = mesh.edge[c->Edge(j)].get();
        Matrix temp = Matrix::Zero();
        Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1];
        for (int k = 0; k <= kOrder; ++k){
          dp[k] = Pow(distance, 2*k-1)/ Pow(Factorial(k), 2);
        }
        e->Integrate([&](const Node& node) {
          return GetMatAt(node.data(), *c, *c, normal, dp);
        }, &temp);
        A_inv[i] += temp;
        e->Integrate([&](const Node& node) {
          return GetVecAt(*c, node.data(), distance);
        }, &(block[offset[i]+j].b_sub));
        // cout << block[offset[i]+j].b_sub.transpose() << endl << endl;
      }
      A_inv[i] = A_inv[i].inverse();
      // cout << A_inv[i] << endl << endl;
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
  vector<int> offset;
  vector<Matrix> A_inv;
  vector<Matrix> B_mat;
  vector<Column> b_col;
  vector<VrBlock> block;

  // functions
  VrApproach() = default;
  void BuildBlock(const MeshType& mesh, const BndCondsType& bndcond);

 private:
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
  static Column GetVecAt(const Cell& a, const Real* coord, Real distance) {
    return a.Functions(coord) / distance;
  }
};

#endif // INCLUDE_VRAPPROACH_HPP_

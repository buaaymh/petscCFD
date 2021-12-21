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
  Array         priVar;
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
  void CalculateBmats(const EdgeGroups& edgeGroups) {
    for (auto& [type, group] : edgeGroups) {
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
  void CalculateAinvs(const Mesh<kOrder>& mesh, const EdgeGroups& edgeGroups) {
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      auto c = mesh.cell[i].get();
      Matrix a_mat = Matrix::Zero();
      for (int j = mesh.offset[i]; j < mesh.offset[i+1]; ++j) {
        auto e = mesh.edge[mesh.edge_csr[j]].get();
        Matrix temp = Matrix::Zero();
        Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1];
        if (mesh.adjc_csr[j] >= 0) { // Interiod, periodic boundary
          a_mat += edgeGroups.at(0)->CalculateMat(e, c, c);
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
        auto e = mesh.edge[mesh.edge_csr[j]].get();
        if (e->right == c) {
          block[j].C_mat = A_inv[i] * B_mat[e->I()];
        } else {
          block[j].C_mat = A_inv[i] * B_mat[e->I()].transpose();
        }
      }
    }
  }
  void UpdatePriVar(const Array& conVar) {
    for (int i = 0; i < conVar.size(); i += nEqual) {
      Physics::ConsToPrim(conVar.data()+i, priVar.data()+i);
    }
  }
  void UpdateBcols(const Mesh<kOrder>& mesh) {
    const auto& offset = mesh.offset;
    const auto& adjc_csr = mesh.adjc_csr;
    for (int i = 0; i < b_col.size(); ++i) {
      b_col[i].setZero();
      for (int j = offset[i]; j < offset[i+1]; ++j) {
        if (adjc_csr[j] >= 0) {
          auto priVar_d = priVar.col(adjc_csr[j]) - priVar.col(i);
          b_col[i] += block[j].b_sub * priVar_d.transpose();
        }
      }
      b_col[i] = A_inv[i] * b_col[i];
    }
  }
  void UpdateCoefs(const Mesh<kOrder>& mesh) {
    const auto& offset = mesh.offset;
    const auto& adjc_csr = mesh.adjc_csr;
    for (int k = 0; k < 7; ++k) {
      for (int i = 0; i < b_col.size(); ++i) {
        EqualCol temp = EqualCol::Zero();
        for (int j = offset[i]; j < offset[i+1]; ++j) {
          if (adjc_csr[j] >= 0) {
            temp += block[j].C_mat * coefs.block<nCoef, nEqual>(0, adjc_csr[j]*nEqual);
          }
        }
        temp += b_col[i];
        coefs.block<nCoef, nEqual>(0, i*nEqual) *= -0.3;
        coefs.block<nCoef, nEqual>(0, i*nEqual) += temp * 1.3;
      }
      if (k % 2 == 0) {
        PetscSFBcastBegin(sfCoef, MPIU_REAL, coefs.data(), coefs.data(), MPI_REPLACE);
        PetscSFBcastEnd(sfCoef, MPIU_REAL, coefs.data(), coefs.data(), MPI_REPLACE);
      }
    }
  }
  void InitLimiter(const Mesh<kOrder>& mesh, const EdgeGroups& edgeGroups);
  void Limiter(const Mesh<kOrder>& mesh);
  void PositivePreserve(const Mesh<kOrder>& mesh);
 
 private:
  vector<Real> dr_sum;
  vector<Real> r_max;
  static Column GetVecAt(const Cell<kOrder>& a, const Real* coord, Real distance) {
    return a.Functions(coord) / distance;
  }
  void LimitTroubleCell(const Mesh<kOrder>& mesh, const Cell<kOrder>* cell);
};

template <int kOrder, class Physics>
void SpaceDiscr<kOrder, Physics>::InitLimiter(const Mesh<kOrder>& mesh,
                                              const EdgeGroups& edgeGroups)
{
  r_max.assign(mesh.CountCells(), EPS);
  dr_sum.assign(mesh.CountCells(), 0.0);
  // For each interior edge
  for (const auto& e : edgeGroups.at(0)->edge) {
    auto cell_l = e->left, cell_r = e->right;
    int i = cell_l->I(), j = cell_r->I();
    // max rho of both sides
    Real rl = priVar(0,i), rr = priVar(0,j);
    Real rMax = max(rl, rr);
    r_max[i] = max(rMax, r_max[i]), r_max[j] = max(rMax, r_max[j]);
    const auto& coef_l = coefs.col(i*nEqual), coef_r = coefs.col(j*nEqual);
    Real dr = Abs((rl + cell_l->Functions(e->Center().data()).dot(coef_l)) -
                  (rr + cell_r->Functions(e->Center().data()).dot(coef_r)));
    dr_sum[i] += dr;
    dr_sum[j] += dr;
  }
  // For each periodic edge
  if (edgeGroups.find(1) != edgeGroups.end()) {
    for (const auto& e : edgeGroups.at(1)->edge) {
      auto cell_l = e->left, cell_r = e->right;
      int i = cell_l->I(), j = cell_r->I();
      // max rho of both sides
      Real rl = priVar(0,i), rr = priVar(0,j);
      Real rMax = max(rl, rr);
      r_max[i] = max(rMax, r_max[i]), r_max[j] = max(rMax, r_max[j]);
      const auto& coef_l = coefs.col(i*nEqual), coef_r = coefs.col(j*nEqual);
      Real dr = Abs((rl + cell_l->Functions(e->Center().data()).dot(coef_l)) -
                    (rr + cell_r->Functions(e->Center().data()).dot(coef_r)));
      dr_sum[i] += dr * 0.5;
      dr_sum[j] += dr * 0.5;
    }
  }
}
template <int kOrder, class Physics>
void SpaceDiscr<kOrder, Physics>::Limiter(const Mesh<kOrder>& mesh)
{
  for (int i = 0; i < mesh.NumLocalCells(); ++i) {
    auto cell = mesh.cell[i].get();
    Real h_k = Pow(cell->Measure(), (kOrder+1)*0.25);
    Real indicator = dr_sum[i] / (cell->nCorner() * h_k * r_max[i]);
    if (indicator > 1) LimitTroubleCell(mesh, cell);
  }
  PetscSFBcastBegin(sfCoef, MPIU_REAL, coefs.data(), coefs.data(), MPI_REPLACE);
  PetscSFBcastEnd(sfCoef, MPIU_REAL, coefs.data(), coefs.data(), MPI_REPLACE);
}
template <int kOrder, class Physics>
void SpaceDiscr<kOrder, Physics>::LimitTroubleCell(const Mesh<kOrder>& mesh,
                                                   const Cell<kOrder>* cell)
{
  const auto& offset = mesh.offset;
  const auto& edge_csr = mesh.edge_csr;
  const auto& adjc_csr = mesh.adjc_csr;
  int i = cell->I();
  const auto& ruvp_l = priVar.col(i);
  Eigen::Array<Real,nEqual,1> d1_max = ruvp_l;
  Eigen::Array<Real,nEqual,1> d1_min = ruvp_l;
  for (int ec = offset[i]; ec < offset[i+1]; ++ec) {
    int j = adjc_csr[ec];
    if (j >= 0) {
      const auto& ruvp_r = priVar.col(j);
      d1_max = d1_max.max(ruvp_r.array());
      d1_min = d1_min.min(ruvp_r.array());
    }
  }
  d1_max -= ruvp_l.array(); d1_min -= ruvp_l.array();
  const auto& coef_l = coefs.block<2, nEqual>(0, i*nEqual);
  Real eps2 = 1 * Pow(cell->Measure(), 1.5);
  Eigen::Matrix<Real,nEqual,1> factor = Eigen::Matrix<Real,nEqual,1>::Ones();
  for (int ec = offset[i]; ec < offset[i+1]; ++ec) {
    auto d2 = cell->Functions(mesh.edge[edge_csr[ec]]->Center().data()).head(2).transpose() * coef_l;
    for (int m = 0; m < nEqual; ++m) {
      if (d2(m) > EPS) {
        Real venkat = Venkat(d1_max(m), d2(m), eps2);
        factor(m) = min(factor(m), venkat);
      } else if (d2(m) < -EPS) {
        Real venkat = Venkat(d1_min(m), d2(m), eps2);
        factor(m) = min(factor(m), venkat);
      }
    }
  }
  for (int m = 0; m < nEqual; ++m) {
    coefs.col(i*nEqual+m).head<2>() *= factor(m);
    coefs.col(i*nEqual+m).tail<nCoef-2>().setZero();
  }
}

template <int kOrder, class Physics>
void SpaceDiscr<kOrder, Physics>::PositivePreserve(const Mesh<kOrder>& mesh)
{
  for (int i = 0; i < mesh.CountCells(); ++i) {
    auto cell = mesh.cell[i].get();
    Node dr{1.0/cell->DxInv(), 1.0/cell->DyInv()};
    dr += cell->Center();
    auto functions = cell->Functions(dr.data());
    auto functions2 = cell->Functions(cell->Center().data());
    // Limit density
    Real u = priVar(0,i);
    Real du = functions.dot(coefs.col(i*nEqual).array().abs().matrix());
    Real du2 = functions2.dot(coefs.col(i*nEqual).array().abs().matrix());
    if (u - du < 0 || u + du2 < 0) {
      Real factor = u / du;
      coefs.col(i*nEqual) *= 0.0;
      coefs.col(i*nEqual+1) *= 0.0;
      coefs.col(i*nEqual+2) *= 0.0;
      coefs.col(i*nEqual+3) *= 0.0;
    }
    // Limit pressure
    u = priVar(3,i);
    du = functions.dot(coefs.col(i*nEqual+3).array().abs().matrix());
    du2 = functions2.dot(coefs.col(i*nEqual+3).array().abs().matrix());
    if (u - du < 0 || u + du2 < 0) {
      Real factor = u / du;
      coefs.col(i*nEqual) *= 0.0;
      coefs.col(i*nEqual+1) *= 0.0;
      coefs.col(i*nEqual+2) *= 0.0;
      coefs.col(i*nEqual+3) *= 0.0;
    }
  }
}

}  // cfd

#endif // INCLUDE_SPACEDISCR_HPP_

/// @file solver.hpp
///
/// Definition of the class related to the flow solver.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_SOLVER_HPP_
#define INCLUDE_SOLVER_HPP_

#include "defs.hpp"
#include "physModel.hpp"
#include "bndConds.hpp"
#include "vrApproach.hpp"
#include "timeDiscr.hpp"
#include "geometry/mesh.hpp"

#include <string>
#include "data/path.hpp"  // defines TEST_DATA_DIR

using namespace std;

template<int kOrder, class Physics>
class Solver
{
 public:
  using MeshType =  Mesh<kOrder>;
  Physics                            physics;     /**< physical model */
  MeshType                           mesh;        /**< element and geometry */
  BndConds<kOrder, Physics>          bndConds;    /**< boundary conditions */
  VrApproach<kOrder, Physics>        vrApproach;  /**< variational reconstruction */
  RK3TS                              ts;

  // functions
  Solver() = default;
  void SetupDataLayout() {
    mesh.SetDataLayout(physics);
    vrApproach.SetCoefLayout(mesh.dm);
  }

  void SetBoundaryConditions(const void* ctx) {
    const BC* bc = (const BC*) ctx;
    bndConds.InitializeBndConds(mesh.dm, bc);
    bndConds.ClassifyEdges(mesh);
    bndConds.PreProcess();
    mesh.UpdateCellNeighbs(bndConds.types);
  }
  void InitializeDS() {
    vrApproach.AllocatorMats(mesh);
    for (auto& [type, bd] : bndConds.bdGroup) {
      bd->CalculateBmats(vrApproach);
    }
    vrApproach.CalculateAinvs(mesh);
    vrApproach.CalculateBlockC(mesh);
  }
  void InitializeTS() {
    ts.SetDM(mesh.dm);
    ts.SetSolverContext(this);
    ts.SetTimeEndAndSetpNum(1.0, 4);
  }
  void InitializeSolution(function<void(int dim, const Real*, int, Real*)>func) {
    Vec         X;
    int         nCell = mesh.CountCells(), Nf = Physics::nEqual;
    Real*       solution;
    Real*       coord;

    DMGetLocalVector(mesh.dm, &X);
    VecGetArray(X, &solution);
    for (int i = 0; i < nCell; ++i) {
      PetscPrintf(PETSC_COMM_SELF, "(%f, %f)\n", mesh.cell[i]->Center()(0), mesh.cell[i]->Center()(1));
      func(2, mesh.cell[i]->Center().data(), Nf, solution+i*Nf);
    }
    VecRestoreArray(X, &solution);
    DMRestoreLocalVector(mesh.dm, &X);
  }
  void CalculateScalar() {
    Vec X;

    ts.SetRHSFunction(RHSFunctionScalar);
    DMGetLocalVector(mesh.dm, &X);
    VecView(X, PETSC_VIEWER_STDOUT_SELF);
    ts.Solver(X);
    VecView(X, PETSC_VIEWER_STDOUT_SELF);
    DMRestoreLocalVector(mesh.dm, &X);
  }


 private:
  static constexpr auto RHSFunctionScalar = [](Real t, Vec U, Vec RHS, void *ctx) {
    using Flux = Eigen::Matrix<Real, Physics::nEqual, 1>;
    Solver<kOrder, Physics>*      solver = (Solver<kOrder, Physics>*) ctx;
    const Real*                   cv;
    Real*                         rhs;
    const int                     nEqual = Physics::nEqual;
    Flux                          fc;

    VecZeroEntries(RHS);
    VecGetArrayRead(U, &cv);
    VecGetArray(RHS, &rhs);

    // for (auto& e : solver->mesh.interior) {
    //   fc = Flux::Zero();
    //   auto cell_l = e->left;
    //   auto cell_r = e->right;
    //   int i = cell_l->I(), j = cell_r->I();
    //   Real normal[2] = {e->Nx(), e->Ny()};
    //   e->Integrate([&](const Node& p) {
    //     Real* u_l = cv[nEqual*i];
    //     Real* u_r = cv[nEqual*j];
    //     return Physics::Riemann(2, p.data(), u_l, u_r, normal);
    //   }, &fc);
    //   for (int k = 0; k < nEqual; ++k) {
    //     rhs[nEqual*i+k] += fc(k);
    //     rhs[nEqual*j+k] -= fc(k);
    //   }
    // }
    // for (auto& [type, bd] : solver->bndConds.bdGroup) {
    //   bd.UpdateRHS(cv, rhs);
    // }

    VecRestoreArrayRead(U, &cv);
    VecRestoreArray(RHS, &rhs);
  };
};

#endif // INCLUDE_SOLVER_HPP_

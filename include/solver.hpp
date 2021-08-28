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

struct User {
  int           order = 2;
  string        filename = "medium.msh";
  string        model = "box";
  Real          cfl = 1.0, tEnd = 1.0, output_interval = 1.0;

  static constexpr auto InitFunc = [](int dim, const Real* coord, int Nf, Real* u) {
    for (int i = 0; i < Nf; i++) {
      u[i] = Sin(2 * coord[0] * PI);
    }
  };
};

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
    for (auto& [type, bd] : bndConds.bdGroup) { bd->PreProcess(); }
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
  void InitializeTS(const void* ctx) {
    const User* user = (const User*) ctx;
    const string dir{OUTPUT_DIR};

    ts.SetSolverContext(this);
    ts.SetTimeEndAndSetpNum(1.0, 20);
    ts.SetOutputInterval(1);
    ts.SetOutputDirModelName(dir + user->model);
  }
  void InitializeSolution(function<void(int dim, const Real*, int, Real*)>func) {
    Vec         X;
    int         nCell = mesh.CountCells(), Nf = Physics::nEqual;
    Real*       solution;

    DMGetLocalVector(mesh.dm, &X);
    VecGetArray(X, &solution);
    for (int i = 0; i < nCell; ++i) {
      func(2, mesh.cell[i]->Center().data(), Nf, solution+i*Nf);
    }
    VecRestoreArray(X, &solution);
    DMRestoreLocalVector(mesh.dm, &X);
  }
  void CalculateScalar() {
    Vec X;
    ts.SetMonitor(OutputScalar);
    ts.SetRHSFunction(RHSFunctionScalar);
    DMGetLocalVector(mesh.dm, &X);
    ts.Solver(mesh.dm, X);
    DMRestoreLocalVector(mesh.dm, &X);
  }

 private:
  static constexpr auto OutputScalar = [](DM dm, Vec U_local,
      const char* filename, PetscViewer viewer) {
    Vec U_global;
    PetscViewerFileSetName(viewer, filename);
    DMGetGlobalVector(dm, &U_global);
    PetscObjectSetName((PetscObject) U_global, "");
    DMLocalToGlobal(dm, U_local, INSERT_VALUES, U_global);
    VecView(U_global, viewer);
    DMRestoreGlobalVector(dm, &U_global);
  };
  static constexpr auto RHSFunctionScalar = [](Real t, Vec U, Vec RHS, void *ctx) {
    using Flux = Eigen::Matrix<Real, Physics::nEqual, 1>;
    Solver<kOrder, Physics>*      solver = static_cast<Solver<kOrder, Physics>*>(ctx);
    const Real*                   cv;
    Real*                         rhs;
    const int                     nEqual = Physics::nEqual;

    VecZeroEntries(RHS);
    VecGetArrayRead(U, &cv);
    VecGetArray(RHS, &rhs);
    for (auto& [type, bd] : solver->bndConds.bdGroup) {
      bd->UpdateRHS(cv, rhs, solver);
    }
    for (auto& c : solver->mesh.cell) {
      int   start = c->I() * nEqual;
      Real  vol = c->Measure();
      for (int i = 0; i < nEqual; ++i) { rhs[start+i] /= vol; }
    }
    PetscSFScatterBegin(solver->mesh.sf, MPIU_REAL, rhs, rhs);
    PetscSFScatterEnd(solver->mesh.sf, MPIU_REAL, rhs, rhs);
    VecRestoreArrayRead(U, &cv);
    VecRestoreArray(RHS, &rhs);
  };
};

#endif // INCLUDE_SOLVER_HPP_

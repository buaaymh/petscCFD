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
#include <iostream>
#include <string>
#include "data/path.hpp"  // defines TEST_DATA_DIR

using namespace std;
using namespace Eigen;

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
  using ConVar = Matrix<Real, Physics::nEqual, Dynamic>;
  Physics                            physics;     /**< physical model */
  MeshType                           mesh;        /**< element and geometry */
  BndConds<kOrder, Physics>          bndConds;    /**< boundary conditions */
  VrApproach<kOrder, Physics>        vrApproach;  /**< variational reconstruction */
  RK3TS<kOrder, Physics>             ts;

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
    int nCell = mesh.CountCells(), Nf = Physics::nEqual;
    bndConds.cv.resize(nCell);
    for (int i = 0; i < nCell; ++i) {
      func(2, mesh.cell[i]->Center().data(), Nf, bndConds.cv.data()+i*Nf);
    }
  }
  void CalculateScalar() {
    ts.SetMonitor(OutputScalar);
    ts.SetRHSFunction(RHSFunction);
    ts.Solver(mesh.dm, bndConds.cv);
  }

 private:
  static constexpr auto OutputScalar = [](DM dm, const ConVar& cv, const char* filename,
                                          PetscViewer viewer)
  {
    int nCell; DMPlexGetDepthStratum(dm, 2, nullptr, &nCell);
    Vec pv_local, pv_global;
    PetscViewerFileSetName(viewer, filename);
    DMGetGlobalVector(dm, &pv_global);
    DMGetLocalVector(dm, &pv_local);
    PetscObjectSetName((PetscObject) pv_global, "");
    if (Physics::nEqual == 1) {
      VecPlaceArray(pv_local, cv.data());
      DMLocalToGlobal(dm, pv_local, INSERT_VALUES, pv_global);
      VecResetArray(pv_local);
    }
    DMRestoreLocalVector(dm, &pv_local);
    VecView(pv_global, viewer);
    DMRestoreGlobalVector(dm, &pv_global);
  };
  static constexpr auto RHSFunction = [](Real t, const ConVar& cv, ConVar& rhs, void *ctx) {
    Solver<kOrder, Physics>*  solver = static_cast<Solver<kOrder, Physics>*>(ctx);
    const int                 nEqual = Physics::nEqual;

    rhs.setZero();
    for (auto& [type, bd] : solver->bndConds.bdGroup) { bd->UpdateRHS(cv, rhs, solver); }
    for (auto& c : solver->mesh.cell) {
      Real vol = c->Measure();
      rhs.col(c->I()) /= vol;
    }
    PetscSFScatterBegin(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data());
    PetscSFScatterEnd(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data());
  };
};

#endif // INCLUDE_SOLVER_HPP_

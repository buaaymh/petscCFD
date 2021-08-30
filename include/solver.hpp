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

#include "defs.h"
#include "physModel.hpp"
#include "edgeTools.hpp"
#include "vrApproach.hpp"
#include "timeDiscr.hpp"
#include "geometry/mesh.hpp"
#include "data/path.hpp"  // defines TEST_DATA_DIR

namespace cfd {

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
  EdgeManager<kOrder, Physics>       edgeManager; /**< boundary conditions */
  VrApproach<kOrder, Physics>        vrApproach;  /**< variational reconstruction */
  RK3TS<kOrder, Physics>             ts;

  // functions
  Solver() = default;
  void SetupDataLayout() {
    mesh.SetDataLayout(physics);
    vrApproach.SetCoefLayout(mesh.dm);
  }

  void SetBoundaryConditions(const void* ctx) {
    const BndConds* bc = (const BndConds*) ctx;
    edgeManager.InitializeBndConds(mesh.dm, bc);
    edgeManager.ClassifyEdges(mesh);
    for (auto& [type, bd] : edgeManager.bdGroup) { bd->PreProcess(); }
    mesh.UpdateCellNeighbs(edgeManager.types);
  }
  void InitializeDS() {
    vrApproach.AllocatorMats(mesh);
    for (auto& [type, bd] : edgeManager.bdGroup) {
      bd->CalculateBmats(vrApproach);
    }
    vrApproach.CalculateAinvs(mesh, edgeManager);
    vrApproach.CalculateBlockC(mesh);
  }
  void InitializeTS(const void* ctx) {
    const User* user = static_cast<const User*>(ctx);
    const string dir{OUTPUT_DIR};
    const string output_dir = dir + user->model + "/";
    system(("rm -rf " + output_dir).c_str());
    system(("mkdir -p " + output_dir).c_str());

    ts.SetSolverContext(this);
    ts.SetTimeEndAndSetpNum(1.0, 20);
    ts.SetOutputInterval(1);
    ts.SetOutputDirModelName(output_dir + user->model);
  }
  void InitializeSolution(function<void(int dim, const Real*, int, Real*)>func) {
    int nCell = mesh.CountCells(), Nf = Physics::nEqual;
    edgeManager.cv.resize(nCell);
    for (int i = 0; i < nCell; ++i) {
      func(2, mesh.cell[i]->Center().data(), Nf, edgeManager.cv.data()+i*Nf);
    }
  }
  void Calculate() {
    ts.SetMonitor(Output);
    ts.SetRHSFunction(RHSFunction);
    ts.Solver(mesh.dm, edgeManager.cv);
  }

 private:
  static constexpr auto Output = [](DM dm, const ConVar& cv, const char* filename,
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
    // Update Coefs for each cell
    
    // Solving Riemann flux and update rhs
    for (auto& [type, bd] : solver->edgeManager.bdGroup) { bd->UpdateRHS(cv, rhs, solver); }
    for (auto& c : solver->mesh.cell) {
      Real vol = c->Measure(); rhs.col(c->I()) /= vol;
    }
    PetscSFScatterBegin(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data());
    PetscSFScatterEnd(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data());
  };
};

}   // cfd

#endif // INCLUDE_SOLVER_HPP_

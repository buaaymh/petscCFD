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

struct User;

template<int kOrder, class Physics>
class Solver
{
 public:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr int nEqual = Physics::nEqual;
  using Vr = typename cfd::VrApproach<kOrder, Physics>;
  using MeshType =  typename cfd::Mesh<kOrder>;
  using ConVar = Matrix<Real, nEqual, Dynamic>;
  Physics                            physics;     /**< physical model */
  MeshType                           mesh;        /**< element and geometry */
  EdgeManager<kOrder, Physics>       edgeManager; /**< boundary conditions */
  Vr                                 vrApproach;  /**< variational reconstruction */
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
    ts.SetTimeEndAndSetpNum(user->tEnd, user->n_step);
    ts.SetOutputInterval(user->output_interval);
    ts.SetOutputDirModelName(output_dir + user->model);
  }
  void InitializeSolution(function<void(int dim, const Real*, int, Real*)>func) {
    int nCell = mesh.CountCells(), Nf = nEqual;
    edgeManager.cv.resize(nEqual, nCell);
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
  static constexpr void UpdateCoefs(const ConVar& cv, Solver<kOrder, Physics>* solver) {
    // Update b_col
    solver->vrApproach.UpdateBcols(solver->mesh, cv, solver->edgeManager);
    // Update coefs
    solver->vrApproach.UpdateCoefs(solver->mesh, solver->edgeManager);

  }
  static constexpr auto Output = [](DM dm, const ConVar& cv, const char* filename,
                                          PetscViewer viewer)
  {
    int nCell; DMPlexGetDepthStratum(dm, 2, nullptr, &nCell);
    Vec pv_local, pv_global;
    PetscViewerFileSetName(viewer, filename);
    DMGetGlobalVector(dm, &pv_global);
    DMGetLocalVector(dm, &pv_local);
    PetscObjectSetName((PetscObject) pv_global, "");
    if (nEqual == 1) {
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

    rhs.setZero();
    // Update Coefs for each cell
    UpdateCoefs(cv, solver);
    // Detect touble cell and limite the Coefs

    // Solving Riemann flux and update rhs
    for (auto& [type, bd] : solver->edgeManager.bdGroup) { bd->UpdateRHS(cv, rhs, solver); }
    for (auto& c : solver->mesh.cell) {
      Real vol = c->Measure(); rhs.col(c->I()) /= vol;
    }
    PetscSFBcastBegin(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data(), MPI_REPLACE);
    PetscSFBcastEnd(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data(), MPI_REPLACE);
  };
};

}   // cfd

#endif // INCLUDE_SOLVER_HPP_

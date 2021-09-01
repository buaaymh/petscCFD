/// @file timeDiscr.hpp
///
/// Definition of the class related to temporal discretization.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_TIMEDISCR_HPP_
#define INCLUDE_TIMEDISCR_HPP_

#include "defs.h"

namespace cfd {

class GlobalRungeKutta
{
 public:
  GlobalRungeKutta() = default;
  /* Before Calculation */
  void SetTimeEndAndSetpNum(Real tEnd, int nStep) {
    tEnd_ = tEnd; nStep_ = nStep;
    dt = tEnd_ / nStep_;
  }
  void SetSolverContext(void* ctx) { ctx_ = ctx; }
  void SetOutputDirModelName(const string& dir_model) { dir_model_ = dir_model; }
  void SetOutputInterval(int interval) { interval_ = interval; }
  void SetMonitor(function<void(DM, const Array&, const char*, PetscViewer)>output) {
    Output = output;
  }
  void SetRHSFunction(function<void(Real, const Array&, Array&, void*)> rhs) { Rhs = rhs; }
  void SetComputeInitialCondition(void(*init) (Vec)) { Init = init; }
  /* During Calculation */
  void Solver(DM dm, Array& conVar) {
    PetscViewer     viewer;
    Real            t_current;
    PetscViewerCreate(PetscObjectComm((PetscObject)dm), &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERVTK);
    auto filename = dir_model_ + "."+ std::to_string(0) + ".vtu";
    Output(dm, conVar, filename.data(), viewer);
    for (int i = 1; i <= nStep_; ++i) {
      TimeStepping(conVar);
      t_current += dt;
      if (i % interval_ == 0) {
        auto filename = dir_model_ + "."+ std::to_string(i) + ".vtu";
        Output(dm, conVar, filename.data(), viewer);
      }
      PetscPrintf(PETSC_COMM_WORLD, "Progress: %D/%D at %.2fs\n", i, nStep_, t_current);
    }
    PetscViewerDestroy(&viewer);
  }

 private:
  void TimeStepping(Array& conVar) {
    Array conVar_stage = conVar;
    Array cd_dt(conVar.rows(), conVar.cols());
    /******** Step 1 ********/
    Rhs(dt, conVar_stage, cd_dt, ctx_);
    conVar_stage += cd_dt * dt;
    /******** Step 2 ********/
    Rhs(dt, conVar_stage, cd_dt, ctx_);
    conVar_stage += cd_dt * dt;
    conVar_stage = 0.25 * conVar_stage + 0.75 * conVar;
    /******** Step 3 ********/
    Rhs(dt, conVar_stage, cd_dt, ctx_);
    conVar_stage += cd_dt * dt;
    conVar = 2.0/3 * conVar_stage + 1.0/3 * conVar;
  }

 private:
  int nStep_, interval_;
  void* ctx_;
  Real tEnd_, dt;
  string dir_model_;
  function<void(Vec)> Init;
  function<void(Real, const Array&, Array&, void*)> Rhs;
  function<void(DM, const Array&, const char*, PetscViewer)> Output;
};

}  // cfd

#endif // INCLUDE_TIMEDISCR_HPP_

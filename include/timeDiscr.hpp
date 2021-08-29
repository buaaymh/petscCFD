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

#include "defs.hpp"
#include <string>
#include <functional>

namespace cfd {

using std::function;

template<int kOrder, class Physics>
class RK3TS
{
 public:
  using ConVar = Eigen::Matrix<Real, Physics::nEqual, Dynamic>;
  RK3TS() = default;
  /* Before Calculation */
  void SetTimeEndAndSetpNum(Real tEnd, int nStep) {
    tEnd_ = tEnd; nStep_ = nStep;
    dt = tEnd_ / nStep_;
  }
  void SetSolverContext(void* ctx) { ctx_ = ctx; }
  void SetOutputDirModelName(const string& dir_model) { dir_model_ = dir_model; }
  void SetOutputInterval(int interval) { interval_ = interval; }
  void SetMonitor(function<void(DM, const ConVar&, const char*, PetscViewer)>output) {
    Output = output;
  }
  void SetRHSFunction(function<void(Real, const ConVar&, ConVar&, void*)> rhs) { Rhs = rhs; }
  void SetComputeInitialCondition(void(*init) (Vec)) { Init = init; }
  /* During Calculation */
  void Solver(DM dm, ConVar& cv) {
    PetscViewer     viewer;
    Real            t_current;
    PetscViewerCreate(PetscObjectComm((PetscObject)dm), &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERVTK);
    auto filename = dir_model_ + "."+ to_string(0) + ".vtu";
    Output(dm, cv, filename.data(), viewer);
    for (int i = 1; i <= nStep_; ++i) {
      TimeStepping(cv);
      t_current += dt;
      if (i % interval_ == 0) {
        auto filename = dir_model_ + "."+ to_string(i) + ".vtu";
        Output(dm, cv, filename.data(), viewer);
      }
      PetscPrintf(PETSC_COMM_WORLD, "Progress: %D/%D at %.2fs\n", i, nStep_, t_current);
    }
    PetscViewerDestroy(&viewer);
  }

 private:
  void TimeStepping(ConVar& cv) {
    ConVar cv_stage = cv;
    ConVar cd_dt(cv.cols());
    /******** Step 1 ********/
    Rhs(dt, cv_stage, cd_dt, ctx_);
    cv_stage += cd_dt * dt;
    /******** Step 2 ********/
    Rhs(dt, cv_stage, cd_dt, ctx_);
    cv_stage += cd_dt * dt;
    cv_stage = 0.25 * cv_stage + 0.75 * cv;
    /******** Step 3 ********/
    Rhs(dt, cv_stage, cd_dt, ctx_);
    cv_stage += cd_dt * dt;
    cv = 2.0/3 * cv_stage + 1.0/3 * cv;
  }

 private:
  int nStep_, interval_;
  void* ctx_;
  Real tEnd_, dt;
  string dir_model_;
  function<void(Vec)> Init;
  function<void(Real, const ConVar&, ConVar&, void*)> Rhs;
  function<void(DM, const ConVar&, const char*, PetscViewer)> Output;
};

}  // cfd

#endif // INCLUDE_TIMEDISCR_HPP_

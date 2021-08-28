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

using namespace std;

class RK3TS
{
 public:
  RK3TS() = default;
  /* Before Calculation */
  void SetTimeEndAndSetpNum(Real tEnd, int nStep) {
    tEnd_ = tEnd; nStep_ = nStep;
    dt = tEnd_ / nStep_;
  }
  void SetSolverContext(void* ctx) { ctx_ = ctx; }
  void SetOutputDirModelName(const string& dir_model) { dir_model_ = dir_model; }
  void SetOutputInterval(int interval) { interval_ = interval; }
  void SetMonitor(function<void(DM, Vec, const char*, PetscViewer)>output) {
    Output = output;
  }
  void SetRHSFunction(function<void(Real, Vec, Vec, void*)> rhs) { Rhs = rhs; }
  void SetComputeInitialCondition(void(*init) (Vec)) { Init = init; }
  /* During Calculation */
  void Solver(DM dm, Vec U) {

    PetscViewer viewer;
    PetscViewerCreate(PetscObjectComm((PetscObject)dm), &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERVTK);
    auto filename = dir_model_ + "."+ to_string(0) + ".vtu";
    Output(dm, U, filename.data(), viewer);
    for (int i = 0; i < nStep_; ++i) {
      TimeStepping(U);
      if (i % interval_ == 0) {
        auto filename = dir_model_ + "."+ to_string(i) + ".vtu";
        Output(dm, U, filename.data(), viewer);
      }
      PetscPrintf(PETSC_COMM_WORLD, "Progress: %D / %D\n", i, nStep_);
    }
    PetscViewerDestroy(&viewer);
  }
  /* After Calculation */

 private:
  void TimeStepping(Vec U) {
    Vec     U_old, RHS;

    VecDuplicate(U, &U_old); VecDuplicate(U, &RHS);
    /******** Step 1 ********/
    VecCopy(U, U_old); 
    Rhs(dt, U_old, RHS, ctx_);
    VecAXPY(U_old, dt, RHS);
    /******** Step 2 ********/
    Rhs(dt, U_old, RHS, ctx_);
    VecAXPY(U_old, dt, RHS);
    VecAXPBY(U_old, 0.75, 0.25, U);
    /******** Step 2 ********/
    Rhs(dt, U_old, RHS, ctx_);
    VecAXPY(U_old, dt, RHS);
    VecAXPBY(U, 2.0/3, 1.0/3, U_old);
    /* Free the space */
    VecDestroy(&U_old);
    VecDestroy(&RHS);
  }

 private:
  int nStep_, interval_;
  void* ctx_;
  Real tEnd_, dt;
  string dir_model_;
  function<void(Vec)> Init;
  function<void(Real, Vec, Vec, void*)> Rhs;
  function<void(DM, Vec, const char*, PetscViewer)> Output;
};

#endif // INCLUDE_TIMEDISCR_HPP_

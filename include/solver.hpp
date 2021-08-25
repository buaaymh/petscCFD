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
#include "petscts.h"
#include <petscsf.h>
#include "geometry/mesh.hpp"

#include <iostream>
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
  BndConds<MeshType, Physics>        bndConds;    /**< boundary conditions */
  VrApproach<kOrder, Physics>        vrApproach;  /**< variational reconstruction */

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
    vrApproach.CalculateBmats(mesh.interior, BdCondType::Interior);
    vrApproach.CalculateAinvs(mesh);
    vrApproach.CalculateBlockC(mesh);
  }
  void InitSolution(void* func);


 private:
  void Output(DM dm, const string& filename)
  {
    Vec               U;
    PetscViewer       viewer;

    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERVTK);
    PetscViewerFileSetName(viewer, filename.data());
    DMGetGlobalVector(dm, &U);
    PetscObjectSetName((PetscObject) U, "");
    VecView(U, viewer);
    PetscViewerDestroy(&viewer);
    DMRestoreGlobalVector(dm, &U);
  }
};

#endif // INCLUDE_SOLVER_HPP_

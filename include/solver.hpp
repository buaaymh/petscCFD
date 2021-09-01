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
#include "spaceDiscr.hpp"
#include "timeDiscr.hpp"
#include "geometry/mesh.hpp"
#include "data/path.hpp"  // defines TEST_DATA_DIR

namespace cfd {

struct User;

template<int kOrder, class Physics>
class Solver
{
 public:
  using EdgeGroups = unordered_map<int, unique_ptr<Group<kOrder,Physics>>>;

  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr int nEqual = Physics::nEqual;

  Physics                            physics;     /**< physical model */
  Mesh<kOrder>                       mesh;        /**< element and geometry */
  unordered_map<int, int>            types;
  EdgeGroups                         faceGroups;
  SpaceDiscr<kOrder, Physics>        spaceDiscr;  /**< variational reconstruction */
  GlobalRungeKutta                   timeStepper;

  // functions
  Solver() = default;
  void SetupDataLayout() {
    mesh.SetDataLayout(physics);
    spaceDiscr.SetCoefLayout(mesh.dm);
  }

  void SetBoundaryConditions(const void* ctx) {
    const BndConds* bc = (const BndConds*) ctx;
    InitializeBndConds(bc);
    ClassifyEdges();
    for (auto& [type, bd] : faceGroups) { bd->PreProcess(); }
    mesh.UpdateCellNeighbs(types);
  }
  void InitializeDS() {
    spaceDiscr.AllocatorMats(mesh);
    spaceDiscr.CalculateBmats(faceGroups);
    spaceDiscr.CalculateAinvs(mesh, faceGroups);
    spaceDiscr.CalculateBlockC(mesh);
  }
  void InitializeTS(const void* ctx) {
    const User* user = static_cast<const User*>(ctx);
    const string dir{OUTPUT_DIR};
    const string output_dir = dir + user->model + "/";
    system(("rm -rf " + output_dir).c_str());
    system(("mkdir -p " + output_dir).c_str());

    timeStepper.SetSolverContext(this);
    timeStepper.SetTimeEndAndSetpNum(user->tEnd, user->n_step);
    timeStepper.SetOutputInterval(user->output_interval);
    timeStepper.SetOutputDirModelName(output_dir + user->model);
  }
  void InitializeSolution(function<void(int dim, const Real*, int, Real*)>func) {
    int nCell = mesh.CountCells(), Nf = nEqual;
    spaceDiscr.conVar.resize(nEqual, nCell);
    for (int i = 0; i < nCell; ++i) {
      func(2, mesh.cell[i]->Center().data(), Nf, spaceDiscr.conVar.data()+i*Nf);
    }
  }
  void Calculate() {
    timeStepper.SetMonitor(Output);
    timeStepper.SetRHSFunction(RHSFunction);
    timeStepper.Solver(mesh.dm, spaceDiscr.conVar);
  }

 private:
  void InitializeBndConds(const BndConds* bc) {
    IS                bdTypeIS;
    DMLabel           label;
    const int         *types;
    int               numTypes;
    /* Interior edge group initialization */
    faceGroups[0] = make_unique<Interior<kOrder,Physics>>();
    /* Boundary edge group initialization */
    DMGetLabel(mesh.dm, "Face Sets", &label);
    DMGetLabelIdIS(mesh.dm, "Face Sets", &bdTypeIS);
    ISGetLocalSize(bdTypeIS, &numTypes);
    ISGetIndices(bdTypeIS, &types);
    for (int i = 0; i < numTypes; ++i) {
      switch (BdCondType(types[i]))
      {
      case BdCondType::Periodic:
        faceGroups[types[i]] = make_unique<Periodic<kOrder,Physics>>(bc->lower, bc->upper);
        break;
      case BdCondType::InFlow:
        faceGroups[types[i]] = make_unique<InFlow<kOrder,Physics>>(bc->inflow);
        break;
      case BdCondType::OutFlow:
        faceGroups[types[i]] = make_unique<OutFlow<kOrder,Physics>>();
        break;
      case BdCondType::FarField:
        faceGroups[types[i]] = make_unique<FarField<kOrder,Physics>>(bc->refVal, 1);
        break;
      case BdCondType::InviscWall:
        faceGroups[types[i]] = make_unique<InviscWall<kOrder,Physics>>();
        break;
      default:
        break;
      }
    }
  }
  void ClassifyEdges() {
    int  eStart, eEnd, type;

    DMPlexGetDepthStratum(mesh.dm, 1, &eStart, &eEnd); /* edges */
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      for(int j = mesh.offset[i]; j < mesh.offset[i+1]; ++j) {
        int e = mesh.edge_csr[j];
        if (mesh.edge[e]->right == nullptr) {
          DMGetLabelValue(mesh.dm, "Face Sets", e+eStart, &type);
          types.emplace(e, type);
          faceGroups[type]->edge.insert(mesh.edge[e].get());
        } else { // interior edges
          faceGroups[0]->edge.insert(mesh.edge[e].get());
        }
      }
    }
  }
  static constexpr void UpdateCoefs(const Array& conVar, Solver<kOrder, Physics>* solver) {
    // Update b_col
    solver->spaceDiscr.UpdateBcols(solver->mesh, conVar);
    // Update coefs
    solver->spaceDiscr.UpdateCoefs(solver->mesh);
    // Detect trouble cells
    // solver->limiter.LimitCoefs(solver->mesh, solver->edgeManager, solver->spaceDiscr);
  }
  static constexpr auto Output = [](DM dm, const Array& conVar, const char* filename,
                                          PetscViewer viewer)
  {
    int nCell; DMPlexGetDepthStratum(dm, 2, nullptr, &nCell);
    Vec pv_local, pv_global;
    PetscViewerFileSetName(viewer, filename);
    DMGetGlobalVector(dm, &pv_global);
    DMGetLocalVector(dm, &pv_local);
    PetscObjectSetName((PetscObject) pv_global, "");
    if (nEqual == 1) {
      VecPlaceArray(pv_local, conVar.data());
      DMLocalToGlobal(dm, pv_local, INSERT_VALUES, pv_global);
      VecResetArray(pv_local);
    }
    DMRestoreLocalVector(dm, &pv_local);
    VecView(pv_global, viewer);
    DMRestoreGlobalVector(dm, &pv_global);
  };
  static constexpr auto RHSFunction = [](Real t, const Array& conVar, Array& rhs, void *ctx) {
    Solver<kOrder, Physics>*  solver = static_cast<Solver<kOrder, Physics>*>(ctx);

    rhs.setZero();
    // Update Coefs for each cell
    UpdateCoefs(conVar, solver);
    // Solving Riemann flux and update rhs
    for (auto& [type, bd] : solver->faceGroups) { bd->UpdateRHS(conVar, solver->spaceDiscr.coefs, rhs); }
    for (auto& c : solver->mesh.cell) {
      Real vol = c->Measure(); rhs.col(c->I()) /= vol;
    }
    PetscSFBcastBegin(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data(), MPI_REPLACE);
    PetscSFBcastEnd(solver->mesh.sf, MPIU_REAL, rhs.data(), rhs.data(), MPI_REPLACE);
  };
};

}   // cfd

#endif // INCLUDE_SOLVER_HPP_

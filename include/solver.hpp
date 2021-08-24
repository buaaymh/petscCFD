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

struct User
{
  int                         outputInterval=1;
  Real                        dt = 1.0, tEnd = 1.0;
  string                      meshfile = "box.msh";
  string                      resultfile = "box.vtu";

  // Periodic Boundary
  Real                        lower[2] = {0.0, 0.0}, upper[2] = {1.0, 1.0};
  // FarField Boundary
  Real                        refVal[1] = {1.0};
  // InFlow Boundary
  static void InFlow(Real t, const Real* coord, Real* Vals) {}
};

template <int kOrder>
class Triangle;
template <int kOrder>
class Quadrangle;

template<int kOrder, class Physics>
class Solver;
template<int kOrder, class Physics>
class Solver
{
 public:
  using MeshType =  Mesh<2, kOrder>;
  static constexpr int               nCoef = (kOrder+1) * (kOrder+2) / 2 - 1;
  TS                                 ts;         /**< temporal discretization */
  MeshType                           mesh;       /**< element and geometry */
  BndConds<MeshType, Physics>        bndConds;   /**< boundary conditions */
  VrApproach<kOrder, Physics>        vrApproach;
  DM                                 dm, dmCoef;
  PetscSF                            sf, sfCoef;
  int                                iter = 0;
  PetscMPIInt                        rank;


  // functions
  Solver(){ MPI_Comm_rank(PETSC_COMM_WORLD, &rank); };
  void ReadMesh(const string filename) {
    DM                   dmDist;
    DM                   dmOverlap;
    PetscPartitioner     part;

    // Read mesh file
    DMPlexCreateFromFile(PETSC_COMM_WORLD, filename.data(), PETSC_TRUE, &dm);
    // Partition mesh
    DMPlexGetPartitioner(dm, &part);
    PetscPartitionerSetType(part,PETSCPARTITIONERPARMETIS);
    DMSetBasicAdjacency(dm, PETSC_TRUE, PETSC_FALSE);
    DMPlexDistribute(dm, 0, nullptr, &dmDist);
    if (dmDist) { DMDestroy(&dm); dm = dmDist; }
    // Record local cell number
    mesh.CountLocalCells(dm);
    // Add overlapped cell
    DMPlexDistributeOverlap(dm, 1, nullptr, &dmOverlap);
    if (dmOverlap) { DMDestroy(&dm); dm = dmOverlap; }
    PetscObjectSetName((PetscObject) dm, "Mesh");
    DMSetFromOptions(dm);
    SetupCoefLayout();
  }
  void SetupDataLayout()
  {
    PetscSection    section;
    int             numFields = Physics::nEqual;
    int             numComp[numFields], numDof[numFields*(DIM+1)];
    int             nroots, nleaves;

    /* Create a PetscSection with this data layout */
    for (int i = 0; i < numFields*(DIM+1); ++i) numDof[i] = 0;
    int i = 0; auto field_desc = Physics::CreateFieldDiscription();
    for (const auto& [name, dof] : field_desc)
    { numDof[i*(DIM+1)+DIM] = dof; numComp[i++] = dof; }
    DMSetNumFields(dm, numFields);
    DMPlexCreateSection(dm, NULL, numComp, numDof, 0, NULL, NULL, NULL, NULL, &section);
    /* Name the Field variables */
    i = 0;
    for (const auto& [name, dof] : field_desc)
      PetscSectionSetFieldName(section, i++, name.c_str());
    PetscSectionSetUp(section);
    /* Tell the DM to use this data layout */
    DMSetLocalSection(dm, section);
    PetscSectionDestroy(&section);
    /* Build the ghosted start forest for data */
    DMGetSectionSF(dm, &sf);
    PetscSFGetGraph(sf, &nroots, &nleaves, nullptr, nullptr);
    int selected[nleaves-nroots];
    for (int i = nroots; i < nleaves; ++i) { selected[i-nroots] = i; }
    PetscSFCreateEmbeddedLeafSF(sf, nleaves-nroots, selected, &sf);
  }
  void ConstructMesh()
  {
    int               cStart, cEnd, eStart, eEnd, vStart, vEnd;
    DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd); /* vertices */
    DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd); /* edges */
    DMPlexGetDepthStratum(dm, 2, &cStart, &cEnd); /* cells */
    mesh.node.reserve(vEnd - vStart);
    mesh.edge.reserve(eEnd - eStart);
    mesh.cell.reserve(cEnd - cStart);
    /* Construct Nodes */
    PetscSection      coordSection;
    Vec               coordinates;
    const Real*       coord;
    DM                dmCoord;
    DMGetCoordinateSection(dm, &coordSection);
    DMGetCoordinatesLocal(dm, &coordinates);
    DMGetCoordinateDM(dm, &dmCoord);
    VecGetArrayRead(coordinates, &coord);
    for (int v = vStart; v < vEnd; ++v) {
      Real*           xy;
      DMPlexPointLocalRead(dmCoord, v, coord, &xy);
      mesh.node.emplace_back(xy[0], xy[1]);
    }
    VecRestoreArrayRead(coordinates, &coord);
    /* Construct Edges */
    for (int e = eStart; e < eEnd; ++e) {
      const int* vPoints;
      DMPlexGetCone(dm, e, &vPoints);
      auto edge_ptr = make_unique<Edge<kOrder>>(e-eStart,
                                                mesh.node.at(vPoints[0]-vStart),
                                                mesh.node.at(vPoints[1]-vStart));
      mesh.edge.emplace_back(move(edge_ptr));
    }
    /* Construct Cells */
    for (int c = cStart; c < cEnd; ++c) {
      int size;
      const int *ePoints, *orientations;
      DMPlexGetConeSize(dm, c, &size); int* corner = new int[size];
      DMPlexGetCone(dm, c, &ePoints);
      DMPlexGetConeOrientation(dm, c, &orientations);
      for (int i = 0; i < size; ++i) {
        const int *vPoints;
        DMPlexGetCone(dm, ePoints[i], &vPoints);
        if (orientations[i] < 0) { corner[i] = vPoints[1]; }
        else { corner[i] = vPoints[0]; }
      }
      if (size == 3) {
        auto triangle_ptr = make_unique<Triangle<kOrder>>(c, mesh.node.at(corner[0]-vStart),
                                                             mesh.node.at(corner[1]-vStart),
                                                             mesh.node.at(corner[2]-vStart));
        triangle_ptr->edge = {ePoints[0]-eStart, ePoints[1]-eStart, ePoints[2]-eStart};
        mesh.cell.emplace_back(move(triangle_ptr));
      }
      if (size == 4) {
        auto quadrangle_ptr = make_unique<Quadrangle<kOrder>>(c, mesh.node.at(corner[0]-vStart),
                                                                 mesh.node.at(corner[1]-vStart),
                                                                 mesh.node.at(corner[2]-vStart),
                                                                 mesh.node.at(corner[3]-vStart));
        quadrangle_ptr->edge = {ePoints[0]-eStart, ePoints[1]-eStart, ePoints[2]-eStart, ePoints[3]-eStart};
        mesh.cell.emplace_back(move(quadrangle_ptr));
      }
      delete[] corner;
      for (int i = 0; i < size; ++i) {
        const int *vPoints;
        DMPlexGetCone(dm, ePoints[i], &vPoints);
        if (orientations[i] < 0) { mesh.edge[ePoints[i]-eStart]->right = mesh.cell[c].get(); }
        else { mesh.edge[ePoints[i]-eStart]->left = mesh.cell[c].get(); }
      }
    }
  }
  void WriteSolution(const void* ctx) {
    const User* user = (const User*) ctx;
    Output(dm, user->resultfile);
  }

  void SetBoundaryConditions(const void* ctx) {
    const User* user = (const User*) ctx;
    BuildBdGroups(user);
    ClassifyBdEdge();
    bndConds.PreProcess();
  }
  void InitializeDS();
  void InitializeTS(const void* ctx)
  {
    const User*       user = (const User*) ctx;
    Vec               U;

    TSCreate(PetscObjectComm((PetscObject)dm), &ts);
    TSSetType(ts,TSSSP); TSSetDM(ts, dm);
    TSSetRHSFunction(ts, nullptr, RHSFunction, this);
    TSSetMaxTime(ts, user->tEnd);
    TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);
    TSSetTimeStep(ts, user->dt); TSSetFromOptions(ts);
    DMGetLocalVector(dm, &U);
    TSSolve(ts,U);
    VecView(U, PETSC_VIEWER_STDOUT_WORLD);
    DMRestoreLocalVector(dm, &U);
  }
  void InitSolution(void* func);


 private:
  void BuildBdGroups(const User* user) {
    IS                bdTypeIS;
    DMLabel           label;
    const int         *types;
    int               numTypes;
    DMGetLabel(dm, "Face Sets", &label);
    DMGetLabelIdIS(dm, "Face Sets", &bdTypeIS);
    ISGetLocalSize(bdTypeIS, &numTypes);
    ISGetIndices(bdTypeIS, &types);
    PetscPrintf(PETSC_COMM_SELF, "**********Processor(%D)**********\n", rank);
    for (int i = 0; i < numTypes; ++i) {
      switch (BdCondType(types[i]))
      {
      case BdCondType::Periodic:
        bndConds.bdGroup[types[i]] = make_unique<PeriodicBd<MeshType,Physics>>(user->lower, user->upper);
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <Periodic>\n", types[i]);
        break;
      case BdCondType::InFlow:
        bndConds.bdGroup[types[i]] = make_unique<InFlowBd<MeshType,Physics>>(user->InFlow);
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <InFlow>\n", types[i]);
        break;
      case BdCondType::OutFlow:
        bndConds.bdGroup[types[i]] = make_unique<OutFlowBd<MeshType,Physics>>();
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <OutFlow>\n", types[i]);
        break;
      case BdCondType::FarField:
        bndConds.bdGroup[types[i]] = make_unique<FarFieldBd<MeshType,Physics>>(user->refVal, 1);
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <FarField>\n", types[i]);
        break;
      default:
        break;
      }
    }
    PetscPrintf(PETSC_COMM_SELF, "******************************\n");
  }
  void ClassifyBdEdge() {
    int               eStart, eEnd;
    DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd); /* edges */
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      for(int j = 0; j < mesh.cell[i]->nCorner(); ++j) {
        int e = mesh.cell[i]->Edge(j);
        if (mesh.edge[e]->left == nullptr || mesh.edge[e]->right == nullptr) {
          int type;
          DMGetLabelValue(dm, "Face Sets", e+eStart, &type);
          bndConds.types.emplace(e, type);
          bndConds.bdGroup[type]->edge.insert(mesh.edge[e].get());
        } else {
          auto dist = (mesh.edge[e]->left->Center() - mesh.edge[e]->right->Center()).norm();
          mesh.edge[e]->SetDist(dist);
          mesh.interior.insert(mesh.edge[e].get());
        }
      }
    }
  }
  void SetupCoefLayout()
  {
    PetscSection    section;
    int             cStart, cEnd;
    int             nroots, nleaves;
    
    DMClone(dm, &dmCoef);
    PetscSectionCreate(PetscObjectComm((PetscObject)dm), &section);
    DMPlexGetHeightStratum(dmCoef, 0, &cStart, &cEnd);
    PetscSectionSetChart(section, cStart, cEnd);
    for (int c = cStart; c < cEnd; ++c) { PetscSectionSetDof(section, c, nCoef); }
    PetscSectionSetUp(section);
    DMSetLocalSection(dmCoef, section);
    PetscSectionDestroy(&section);
    DMGetSectionSF(dmCoef, &sfCoef);
    /* Build the ghosted start forest for data */
    DMGetSectionSF(dmCoef, &sfCoef);
    PetscSFGetGraph(sfCoef, &nroots, &nleaves, nullptr, nullptr);
    int selected[nleaves-nroots];
    for (int i = nroots; i < nleaves; ++i) { selected[i-nroots] = i; }
    PetscSFCreateEmbeddedLeafSF(sfCoef, nleaves-nroots, selected, &sfCoef);
  }
  void Output(DM dm, const string name)
  {
    Vec               U;
    PetscViewer       viewer;
    string            filename = string(OUTPUT_DIR) + name;

    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERVTK);
    PetscViewerFileSetName(viewer, filename.data());
    DMGetGlobalVector(dm, &U);
    PetscObjectSetName((PetscObject) U, "");
    VecView(U, viewer);
    PetscViewerDestroy(&viewer);
    DMRestoreGlobalVector(dm, &U);
  }

  static PetscErrorCode RHSFunction(TS ts, Real dt, Vec U_local, Vec R_local, void *ctx)
  {
    Solver<kOrder, Physics>*     solver = (Solver<kOrder, Physics>*) ctx;

    VecZeroEntries(R_local);
    PetscPrintf(PETSC_COMM_WORLD, "Stepping: %D \n", (solver->iter)++);
    VecScatterBegin(solver->sf, R_local, R_local, INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(solver->sf, R_local, R_local, INSERT_VALUES,SCATTER_FORWARD);
  }
};

#endif // INCLUDE_SOLVER_HPP_

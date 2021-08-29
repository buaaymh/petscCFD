/// @file mesh.hpp
///
/// Definition of the class related to mesh and to metrics.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_MESH_HPP_
#define INCLUDE_MESH_HPP_

#include "defs.hpp"
#include "geometry/element.hpp"
#include <petscsf.h>
#include <memory>
#include <vector>

namespace cfd {

template <int kOrder>
class Mesh
{
 public:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1;
  using EdgeType = Edge<kOrder>;
  using CellType = Cell<kOrder>;
  using GhostType = Ghost<kOrder>;
  
  DM                                 dm;
  PetscSF                            sf;

  std::vector<Node> node;                 /**< node(local+overlap) list */
  std::vector<std::unique_ptr<EdgeType>> edge; /**< edge(local+overlap) list */
  std::vector<std::unique_ptr<CellType>> cell; /**< cell(local+overlap+ghost) list */

  // functions
  Mesh() = default;

  void ReadMeshFile(const std::string& filename) {
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
    DMPlexGetDepthStratum(dm, 2, nullptr, &nCellroot);
    // Add overlapped cell
    DMPlexDistributeOverlap(dm, 1, nullptr, &dmOverlap);
    if (dmOverlap) { DMDestroy(&dm); dm = dmOverlap; }
    PetscObjectSetName((PetscObject) dm, "Mesh");
    DMSetFromOptions(dm);
    ConstructElements();
  }
  template <class Physics>
  void SetDataLayout(const Physics& physics) {
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
  void UpdateCellNeighbs(const std::unordered_map<int, int>& dict) {
    for (int i = 0; i < nCellroot; ++i) {
      for (int j = 0; j < cell[i]->nCorner(); ++j) {
        auto e = edge[cell[i]->Edge(j)].get();
        if (e->left == cell[i].get()) {
          if (e->right) {
            cell[i]->SetAdjc(j, e->right->I());
          } else {
            int type = dict.at(e->I());
            cell[i]->SetAdjc(j, -type);}
        } else {
          cell[i]->SetAdjc(j, e->left->I());
        }
      }
    }
  }
  int NumLocalCells() const { return nCellroot; }
  void Clear() { node.clear(); edge.clear(); cell.clear(); }
  int CountNodes() const { return node.size(); }
  int CountEdges() const { return edge.size(); }
  int CountCells() const { return cell.size(); }
  ~Mesh() = default;

 private:
  Mesh( const Mesh &mesh ) = default;   // override default copy constructor
  Mesh & operator = (const Mesh &mesh); // and assignment operator
  void ConstructElements() {
    int cStart, cEnd, eStart, eEnd, vStart, vEnd;
    DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd); /* vertices */
    DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd); /* edges */
    DMPlexGetDepthStratum(dm, 2, &cStart, &cEnd); /* cells */
    node.reserve(vEnd - vStart);
    edge.reserve(eEnd - eStart);
    cell.reserve(cEnd - cStart);
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
      node.emplace_back(xy[0], xy[1]);
    }
    VecRestoreArrayRead(coordinates, &coord);
    /* Construct Edges */
    for (int e = eStart; e < eEnd; ++e) {
      const int* vPoints;
      DMPlexGetCone(dm, e, &vPoints);
      auto edge_ptr = std::make_unique<EdgeType>(e-eStart,
                                            node.at(vPoints[0]-vStart),
                                            node.at(vPoints[1]-vStart));
      edge.emplace_back(move(edge_ptr));
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
        auto triangle_ptr = std::make_unique<Triangle<kOrder>>(c, node.at(corner[0]-vStart),
                                                             node.at(corner[1]-vStart),
                                                             node.at(corner[2]-vStart));
        triangle_ptr->edge = {ePoints[0]-eStart, ePoints[1]-eStart, ePoints[2]-eStart};
        cell.emplace_back(std::move(triangle_ptr));
      }
      if (size == 4) {
        auto quadrangle_ptr = std::make_unique<Quadrangle<kOrder>>(c, node.at(corner[0]-vStart),
                                                                 node.at(corner[1]-vStart),
                                                                 node.at(corner[2]-vStart),
                                                                 node.at(corner[3]-vStart));
        quadrangle_ptr->edge = {ePoints[0]-eStart, ePoints[1]-eStart, ePoints[2]-eStart, ePoints[3]-eStart};
        cell.emplace_back(std::move(quadrangle_ptr));
      }
      delete[] corner;
      for (int i = 0; i < size; ++i) {
        const int *vPoints;
        DMPlexGetCone(dm, ePoints[i], &vPoints);
        if (orientations[i] < 0) { edge[ePoints[i]-eStart]->right = cell[c].get(); }
        else { edge[ePoints[i]-eStart]->left = cell[c].get(); }
      }
    }
  }
  int nCellroot;
};

}  // cfd

#endif // INCLUDE_MESH_HPP_

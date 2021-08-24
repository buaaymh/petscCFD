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

#include <memory>
#include <set>
#include <petscdmplex.h>

using namespace std;
using Node = Eigen::Matrix<Real,2,1>;

template <int kDim, int kOrder>
class Mesh;

template <int kOrder>
class Mesh<2, kOrder>
{
 public:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1;
  using EdgeType = Edge<kOrder>;
  using CellType = Cell<kOrder>;
  using GhostType = Ghost<kOrder>;

  vector<Node> node;                 /**< node(local+overlap) list */
  vector<unique_ptr<EdgeType>> edge; /**< edge(local+overlap) list */
  vector<unique_ptr<CellType>> cell; /**< cell(local+overlap+ghost) list */

  set<EdgeType*> interior;

  // functions
  Mesh() = default;
  int NumLocalCells() const { return nCellroot; }
  void Clear() { node.clear(); edge.clear(); cell.clear(); }
  void CountLocalCells(const DM& dm) {
    DMPlexGetDepthStratum(dm, 2, nullptr, &nCellroot); /* cells */
  }
  int CountNodes() const { return node.size(); }
  int CountEdges() const { return edge.size(); }
  int CountCells() const { return cell.size(); }
  ~Mesh() = default;

 private:
  Mesh( const Mesh &mesh ) = default;   // override default copy constructor
  Mesh & operator = (const Mesh &mesh); // and assignment operator

  int nCellroot;
};

#endif // INCLUDE_MESH_HPP_

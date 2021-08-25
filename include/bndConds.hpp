/// @file bndConds.hpp
///
/// Definition of the class related to boundary conditions.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_BNDCONDS_HPP_
#define INCLUDE_BNDCONDS_HPP_

#include "defs.hpp"
#include <set>
#include <algorithm>
#include <iostream>
#include "geometry/mesh.hpp"
#include "solver.hpp"

using namespace std;

/// Types of fluid flow
enum class BdCondType { Interior, Periodic, InFlow, OutFlow,
                        FarField, InviscWall, Symmetry };

struct BC {
  // Periodic Boundary
  Real                        lower[2] = {0.0, 0.0}, upper[2] = {1.0, 1.0};
  // FarField Boundary
  Real                        refVal[1] = {1.0};
  // InFlow Boundary
  static void InFlow(Real t, const Real* coord, Real* Vals) {}
};

template<class Mesh, class Physics>
struct Bd {
  using Edge = typename Mesh::EdgeType;
  using Cell = typename Mesh::CellType;
  using Set = typename Mesh::EdgeSet;
  Bd() = default;
  virtual void PreProcess() = 0;
  Set edge;
};
template<class Mesh, class Physics>
struct PeriodicBd : public Bd<Mesh,Physics> {
  using Edge = typename Mesh::EdgeType;
  using Cell = typename Mesh::CellType;
  using Ghost = typename Mesh::GhostType;
  PeriodicBd(const Real* lower, const Real* upper) : Bd<Mesh,Physics>() {
    lower_[0] = lower[0]; lower_[1] = lower[1];
    upper_[0] = upper[0]; upper_[1] = upper[1];
  }
  void PreProcess() override {
    for ( auto& e : Bd<Mesh,Physics>::edge) {
      if (Abs(e->Center()(0)-lower_[0]) < 1e-6) { bdL.emplace_back(e); }
      else if (Abs(e->Center()(1)-lower_[1]) < 1e-6) { bdB.emplace_back(e); }
      else if (Abs(e->Center()(0)-upper_[0]) < 1e-6) { bdR.emplace_back(e); }
      else if (Abs(e->Center()(1)-upper_[1]) < 1e-6) { bdT.emplace_back(e); }
    }
    auto cmpY = [](Edge* a, Edge* b) {
      auto c_a = a->Center()(1), c_b = b->Center()(1);
      return c_a < c_b; };
    auto cmpX = [](Edge* a, Edge* b) {
      auto c_a = a->Center()(0), c_b = b->Center()(0);
      return c_a < c_b; };
    sort(bdL.begin(), bdL.end(), cmpY);
    sort(bdR.begin(), bdR.end(), cmpY);
    sort(bdB.begin(), bdB.end(), cmpX);
    sort(bdT.begin(), bdT.end(), cmpX);
    ghost.reserve(Bd<Mesh,Physics>::edge.size());
    for (int i = 0; i < bdL.size(); ++i) { MatchEdges(bdL[i], bdR[i]); }
    for (int i = 0; i < bdB.size(); ++i) { MatchEdges(bdB[i], bdT[i]); }
  }
  void MatchEdges(Edge* a, Edge* b) {
    auto ab = Node(a->Center() - b->Center());
    if (a->left == nullptr) {
      if (b->right == nullptr) {
        auto b_right = make_unique<Ghost>(*(a->right));
        auto a_left = make_unique<Ghost>(*(b->left));
        a_left->Move(ab); b_right->Move(-ab);
        a->left = a_left.get(); b->right = b_right.get();
        ghost.emplace_back(move(a_left));
        ghost.emplace_back(move(b_right));
      } else {
        auto b_left = make_unique<Ghost>(*(a->right));
        auto a_left = make_unique<Ghost>(*(b->right));
        a_left->Move(ab); b_left->Move(-ab);
        a->left = a_left.get(); b->left = b_left.get();
        ghost.emplace_back(move(a_left));
        ghost.emplace_back(move(b_left));
      }
    } else {
      if (b->right == nullptr) {
        auto b_right = make_unique<Ghost>(*(a->left));
        auto a_right = make_unique<Ghost>(*(b->left));
        a_right->Move(ab); b_right->Move(-ab);
        a->right = a_right.get(); b->right = b_right.get();
        ghost.emplace_back(move(a_right));
        ghost.emplace_back(move(b_right));
      } else {
        auto b_left = make_unique<Ghost>(*(a->left));
        auto a_right = make_unique<Ghost>(*(b->right));
        a_right->Move(ab); b_left->Move(-ab);
        a->right = a_right.get(); b->left = b_left.get();
        ghost.emplace_back(move(a_right));
        ghost.emplace_back(move(b_left));
      }
    }
    auto dist = (a->left->Center() - a->right->Center()).norm();
    a->SetDist(dist); b->SetDist(dist);
  }
  vector<unique_ptr<Cell>> ghost;
  vector<Edge*> bdL;
  vector<Edge*> bdR;
  vector<Edge*> bdB;
  vector<Edge*> bdT;
  Real lower_[2];
  Real upper_[2];
};
template<class Mesh, class Physics>
struct OutFlowBd : public Bd<Mesh,Physics> {
  using Edge = typename Mesh::EdgeType;
  using Cell = typename Mesh::CellType;
  OutFlowBd() : Bd<Mesh,Physics>() {}
  void PreProcess() override {
    for ( auto& e : Bd<Mesh,Physics>::edge) {
      if (e->left == nullptr) { e->left = e->right; }
      else { e->right = e->left; }
      auto dist = (e->left->Center()-e->Center()).norm();
      e->SetDist(dist);
    }
  }
};
template<class Mesh, class Physics>
struct InFlowBd : public Bd<Mesh,Physics> {
  using Value = Eigen::Matrix<Real, Physics::nEqual, 1>;
  using Coefs = Eigen::Matrix<Real, Mesh::nCoef*Physics::nEqual, 1>;
  InFlowBd(void(*func) (Real, const Real*, Real*)) : Bd<Mesh,Physics>(), func_(func) {}
  void PreProcess() override {
    for (auto& e : Bd<Mesh,Physics>::edge) {
      bdVal.emplace(e->I(), Value::Zero());
      bdCoef.emplace(e->I(), Coefs::Zero());
      if (e->left == nullptr) {
        Real dist = (e->right->Center() - e->Center()).norm();
        e->SetDist(dist);
      } else {
        Real dist = (e->left->Center() - e->Center()).norm();
        e->SetDist(dist);
      }
    }
  }
  void(*func_) (Real, const Real*, Real*);
  unordered_map<int, Value> bdVal;
  unordered_map<int, Coefs> bdCoef;
};
template<class Mesh, class Physics>
struct FarFieldBd : public Bd<Mesh,Physics> {
  using Edge = typename Mesh::EdgeType;
  using Cell = typename Mesh::CellType;
  FarFieldBd(const Real* reference, int nVals) : Bd<Mesh,Physics>() {
    refer_ = new Real[nVals];
    for (int i = 0; i < nVals; ++i) { refer_[i] = reference[i]; }
  }
  void PreProcess() override {
    for (auto& e : Bd<Mesh,Physics>::edge) {
      if (e->left == nullptr) {
        Real dist = (e->right->Center() - e->Center()).norm();
        e->SetDist(dist);
      } else {
        Real dist = (e->left->Center() - e->Center()).norm();
        e->SetDist(dist);
      }
    }
  }
  ~FarFieldBd() { delete[] refer_;}
  Real* refer_;
};
template<class Mesh, class Physics>
struct InviscWallBd : public Bd<Mesh,Physics> {
  using Edge = typename Mesh::EdgeType;
  using Cell = typename Mesh::CellType;
  InviscWallBd() : Bd<Mesh,Physics>() {}
  void PreProcess() override {
    for (auto& e : Bd<Mesh,Physics>::edge) {
      if (e->left == nullptr) {
        Real dist = (e->right->Center() - e->Center()).norm();
        e->SetDist(dist);
      } else {
        Real dist = (e->left->Center() - e->Center()).norm();
        e->SetDist(dist);
      }
    }
  }
};
template<class Mesh, class Physics>
struct BndConds {
  BndConds() = default;
  void PreProcess() {
    for (auto& [type, bd] : bdGroup) {
      bd->PreProcess();
    }
  }
  void InitializeBndConds(DM dm, const BC* bc) {
    PetscMPIInt       rank;
    IS                bdTypeIS;
    DMLabel           label;
    const int         *types;
    int               numTypes;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    DMGetLabel(dm, "Face Sets", &label);
    DMGetLabelIdIS(dm, "Face Sets", &bdTypeIS);
    ISGetLocalSize(bdTypeIS, &numTypes);
    ISGetIndices(bdTypeIS, &types);
    PetscPrintf(PETSC_COMM_SELF, "**********Processor(%D)**********\n", rank);
    for (int i = 0; i < numTypes; ++i) {
      switch (BdCondType(types[i]))
      {
      case BdCondType::Periodic:
        bdGroup[types[i]] = make_unique<PeriodicBd<Mesh,Physics>>(bc->lower, bc->upper);
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <Periodic>\n", types[i]);
        break;
      case BdCondType::InFlow:
        bdGroup[types[i]] = make_unique<InFlowBd<Mesh,Physics>>(bc->InFlow);
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <InFlow>\n", types[i]);
        break;
      case BdCondType::OutFlow:
        bdGroup[types[i]] = make_unique<OutFlowBd<Mesh,Physics>>();
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <OutFlow>\n", types[i]);
        break;
      case BdCondType::FarField:
        bdGroup[types[i]] = make_unique<FarFieldBd<Mesh,Physics>>(bc->refVal, 1);
        PetscPrintf(PETSC_COMM_SELF, "Boundary(%D) <FarField>\n", types[i]);
        break;
      default:
        break;
      }
    }
    PetscPrintf(PETSC_COMM_SELF, "******************************\n");
  }
  void ClassifyEdges(Mesh& mesh) {
    int               eStart, eEnd;
    DMPlexGetDepthStratum(mesh.dm, 1, &eStart, &eEnd); /* edges */
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      for(int j = 0; j < mesh.cell[i]->nCorner(); ++j) {
        int e = mesh.cell[i]->Edge(j);
        if (mesh.edge[e]->left == nullptr || mesh.edge[e]->right == nullptr) {
          int type;
          DMGetLabelValue(mesh.dm, "Face Sets", e+eStart, &type);
          types.emplace(e, type);
          bdGroup[type]->edge.insert(mesh.edge[e].get());
        } else {
          auto dist = (mesh.edge[e]->left->Center() - mesh.edge[e]->right->Center()).norm();
          mesh.edge[e]->SetDist(dist);
          mesh.interior.insert(mesh.edge[e].get());
        }
      }
    }
  }
  unordered_map<int, int> types;
  unordered_map<int, unique_ptr<Bd<Mesh,Physics>>> bdGroup;
  ~BndConds() = default;
};
#endif // INCLUDE_BNDCONDS_HPP_

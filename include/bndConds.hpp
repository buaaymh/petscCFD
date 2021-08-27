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
#include "solver.hpp"
#include "vrApproach.hpp"
#include <set>
#include <iostream>

using namespace std;

struct BC {
  // Periodic Boundary
  Real                        lower[2] = {0.0, 0.0}, upper[2] = {1.0, 1.0};
  // FarField Boundary
  Real                        refVal[1] = {1.0};
  // InFlow Boundary
  static void InFlow(Real t, const Real* coord, Real* Vals) {}
};

template<int kOrder, class Physics>
struct EdgeGroup {
  using SolverType = Solver<kOrder, Physics>;
  using Mesh = typename SolverType::MeshType;
  using Edge = typename Mesh::EdgeType;
  struct cmp {bool operator()(Edge* a, Edge* b) const {return a->I() < b->I();}};
  using EdgeSet = set<Edge*, cmp>;
  virtual void PreProcess() = 0;
  virtual void CalculateBmats(VrApproach<kOrder, Physics>& vr) const = 0;
  EdgeSet edge;
};
template<int kOrder, class Physics>
struct Interior : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  // Construction:
  Interior() : EdgeGroup<kOrder,Physics>() {}
  void PreProcess() override {
    for (auto& e : Base::edge) {
      auto dist = (e->left->Center() - e->right->Center()).norm();
      e->SetDist(dist);
    }
  }
  void CalculateBmats(VrApproach<kOrder, Physics>& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
      Real dp[kOrder+1]; Dp<kOrder>::GetDpArrayInterior(distance, dp);
      e->Integrate([&](const Node& node) {
        return vr.GetMatAt(node.data(), *(e->left), *(e->right), normal, dp);
      }, &vr.B_mat[e->I()]);
    }
  }
};
template<int kOrder, class Physics>
struct Periodic : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  using SolverType = Solver<kOrder, Physics>;
  using Mesh = typename SolverType::MeshType;
  using Edge = typename Mesh::EdgeType;
  using Cell = typename Mesh::CellType;
  using Ghost = typename Mesh::GhostType;
  // Construction:
  Periodic(const Real* lower, const Real* upper) : Base() {
    lower_[0] = lower[0]; lower_[1] = lower[1];
    upper_[0] = upper[0]; upper_[1] = upper[1];
  }
  void PreProcess() override {
    for (auto& e : Base::edge) {
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
    ghost.reserve(Base::edge.size());
    for (int i = 0; i < bdL.size(); ++i) { MatchEdges(bdL[i], bdR[i]); }
    for (int i = 0; i < bdB.size(); ++i) { MatchEdges(bdB[i], bdT[i]); }
  }
  void CalculateBmats(VrApproach<kOrder, Physics>& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
      Real dp[kOrder+1]; Dp<kOrder>::GetDpArrayInterior(distance, dp);
      e->Integrate([&](const Node& node) {
        return vr.GetMatAt(node.data(), *(e->left), *(e->right), normal, dp);
      }, &vr.B_mat[e->I()]);
    }
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
  vector<Edge*> bdL, bdR, bdB, bdT;
  Real lower_[2];
  Real upper_[2];
};
template<int kOrder, class Physics>
struct OutFlow : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  // Construction:
  OutFlow() : Base() {}
  void PreProcess() override {
    for ( auto& e : Base::edge) {
      if (e->left == nullptr) { e->left = e->right; }
      else { e->right = e->left; }
      auto dist = (e->left->Center() - e->Center()).norm();
      e->SetDist(dist);
    }
  }
  void CalculateBmats(VrApproach<kOrder, Physics>& vr) const {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
      Real dp[kOrder+1]; Dp<kOrder>::GetDpArrayInterior(distance, dp);
      e->Integrate([&](const Node& node) {
        return vr.GetMatAt(node.data(), *(e->left), *(e->right), normal, dp);
      }, &vr.B_mat[e->I()]);
    }
  }
};
template<int kOrder, class Physics>
struct InFlow : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  using SolverType = Solver<kOrder, Physics>;
  using Mesh = typename SolverType::MeshType;
  using State = Eigen::Matrix<Real, Physics::nEqual, 1>;
  using Coefs = Eigen::Matrix<float, Mesh::nCoef * Physics::nEqual, 1>;
  // Construction:
  InFlow(void(*func) (Real, const Real*, Real*)) : Base(), func_(func) {}
  void PreProcess() override {
    for (auto& e : Base::edge) {
      bdState.emplace(e->I(), State::Zero());
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
  void CalculateBmats(VrApproach<kOrder, Physics>& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1]; Dp<kOrder>::GetDpArrayP0(distance, dp);
        auto c = e->left; if (e->right) c = e->right;
        e->Integrate([&](const Node& node) {
          return vr.GetMatAt(node.data(), *c, *c, normal, dp);
        }, &vr.B_mat[e->I()]);
    }
  }
  void(*func_) (Real, const Real*, Real*);
  unordered_map<int, State> bdState;
  unordered_map<int, Coefs> bdCoef;
};
template<int kOrder, class Physics>
struct FarField : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  FarField(const Real* reference, int nVals) : Base() {
    refer_ = new Real[nVals];
    for (int i = 0; i < nVals; ++i) { refer_[i] = reference[i]; }
  }
  void PreProcess() override {
    for (auto& e : Base::edge) {
      if (e->left == nullptr) {
        Real dist = (e->right->Center() - e->Center()).norm();
        e->SetDist(dist);
      } else {
        Real dist = (e->left->Center() - e->Center()).norm();
        e->SetDist(dist);
      }
    }
  }
  void CalculateBmats(VrApproach<kOrder, Physics>& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1]; Dp<kOrder>::GetDpArrayP0(distance, dp);
        auto c = e->left; if (e->right) c = e->right;
        e->Integrate([&](const Node& node) {
          return vr.GetMatAt(node.data(), *c, *c, normal, dp);
        }, &vr.B_mat[e->I()]);
    }
  }
  ~FarField() { delete[] refer_;}
  Real* refer_;
};
template<int kOrder, class Physics>
struct InviscWall : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  InviscWall() : Base() {}
  void PreProcess() override {
    for (auto& e : Base::edge) {
      if (e->left == nullptr) {
        Real dist = (e->right->Center() - e->Center()).norm();
        e->SetDist(dist);
      } else {
        Real dist = (e->left->Center() - e->Center()).norm();
        e->SetDist(dist);
      }
    }
  }
  void CalculateBmats(VrApproach<kOrder, Physics>& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1]; Dp<kOrder>::GetDpArrayP0(distance, dp);
        auto c = e->left; if (e->right) c = e->right;
        e->Integrate([&](const Node& node) {
          return vr.GetMatAt(node.data(), *c, *c, normal, dp);
        }, &vr.B_mat[e->I()]);
    }
  }
};
template<int kOrder, class Physics>
struct BndConds {
  // Type:
  using SolverType = Solver<kOrder, Physics>;
  using Mesh = typename SolverType::MeshType;
  BndConds() = default;
  void PreProcess() {
    for (auto& [type, bd] : bdGroup) { bd->PreProcess(); }
  }
  void InitializeBndConds(DM dm, const BC* bc) {
    IS                bdTypeIS;
    DMLabel           label;
    const int         *types;
    int               numTypes;
    /* Interior edge group initialization */
    bdGroup[0] = make_unique<Interior<kOrder,Physics>>();
    /* Boundary edge group initialization */
    DMGetLabel(dm, "Face Sets", &label);
    DMGetLabelIdIS(dm, "Face Sets", &bdTypeIS);
    ISGetLocalSize(bdTypeIS, &numTypes);
    ISGetIndices(bdTypeIS, &types);
    for (int i = 0; i < numTypes; ++i) {
      switch (BdCondType(types[i]))
      {
      case BdCondType::Periodic:
        bdGroup[types[i]] = make_unique<Periodic<kOrder,Physics>>(bc->lower, bc->upper);
        break;
      case BdCondType::InFlow:
        bdGroup[types[i]] = make_unique<InFlow<kOrder,Physics>>(bc->InFlow);
        break;
      case BdCondType::OutFlow:
        bdGroup[types[i]] = make_unique<OutFlow<kOrder,Physics>>();
        break;
      case BdCondType::FarField:
        bdGroup[types[i]] = make_unique<FarField<kOrder,Physics>>(bc->refVal, 1);
        break;
      default:
        break;
      }
    }
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
        } else { bdGroup[0]->edge.insert(mesh.edge[e].get()); }
      }
    }
  }
  unordered_map<int, int> types;
  unordered_map<int, unique_ptr<EdgeGroup<kOrder,Physics>>> bdGroup;
  ~BndConds() = default;
};
#endif // INCLUDE_BNDCONDS_HPP_

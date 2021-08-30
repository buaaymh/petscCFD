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

#include "defs.h"

namespace cfd {

template<int kOrder, class Physics>
struct EdgeGroup {
  // Constants:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr int nEqual = Physics::nEqual;
  // Types:
  using Solver = typename cfd::Solver<kOrder, Physics>;
  using Vr = VrApproach<kOrder, Physics>;
  using Mesh = typename cfd::Mesh<kOrder>;
  using Edge = typename Mesh::EdgeType;
  struct cmp {bool operator()(Edge* a, Edge* b) const {return a->I() < b->I();}};
  using EdgeSet = std::set<Edge*, cmp>;
  using ConVar = typename Physics::ConVar;
  using Coefs = typename Eigen::Matrix<Real, Dynamic, Dynamic>;
  using Flux = Matrix<Real, nEqual, 1>;
  // Functions:
  virtual void PreProcess() {
    for (auto& e : edge) {
      auto dist = (e->left->Center() - e->Center()).norm();
      e->SetDist(dist);
    }
  }
  virtual void UpdateVrbCol(const Real* cv, const Mesh& mesh, Vr& vr) const {}
  virtual void CalculateBmats(Vr& vr) const {
    for (auto& e : edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
      Real dp[kOrder+1]; InteriorDp(kOrder, distance, dp);
      e->Integrate([&](const Node& node) {
        return vr.GetMatAt(node.data(), *(e->left), *(e->right), normal, dp);
      }, &vr.B_mat[e->I()]);
    }
  }
  static void InteriorRHS(const Edge* e, const ConVar& cv, const Coefs& coefs, ConVar& rhs) {
    auto cell_l = e->left, cell_r = e->right;
    Flux fc = Flux::Zero();
    const Real normal[2] = {e->Nx(), e->Ny()};
    int i = cell_l->I(), j = cell_r->I();
    e->Integrate([&](const Node& p){
      Flux U_l = cv.col(i) + cell_l->Functions(p.data()).transpose() *
                 coefs.block<nCoef, nEqual>(0, i*nEqual);
      Flux U_r = cv.col(j) + cell_r->Functions(p.data()).transpose() *
                 coefs.block<nCoef, nEqual>(0, j*nEqual);
      return Physics::Riemann(2, p.data(), normal, U_l.data(),
                                                   U_r.data());
    }, &fc);
    rhs.col(cell_l->I()) -= fc;
    rhs.col(cell_r->I()) += fc;
  }
  virtual void UpdateRHS(const ConVar& cv, ConVar& rhs, void* ctx) const {
    Solver*  solver = static_cast<Solver*>(ctx);
    for (auto& e : edge) { InteriorRHS(e, cv, solver->vrApproach.coefs, rhs); }
  }
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
};

template<int kOrder, class Physics>
struct Periodic : public EdgeGroup<kOrder,Physics> {
  // Constants:
  static constexpr int nEqual = Physics::nEqual;
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  using Solver = typename cfd::Solver<kOrder, Physics>;
  using Mesh = typename Solver::MeshType;
  using Edge = typename Mesh::EdgeType;
  using Cell = typename Mesh::CellType;
  using ConVar = typename Physics::ConVar;
  using Flux = Eigen::Matrix<Real, nEqual, 1>;
  // Functions:
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
  virtual void UpdateRHS(const ConVar& cv, ConVar& rhs, void* ctx) const {
    const Solver* solver = static_cast<const Solver*>(ctx);
    for (auto& e : bdL) { Base::InteriorRHS(e, cv, solver->vrApproach.coefs, rhs); }
    for (auto& e : bdB) { Base::InteriorRHS(e, cv, solver->vrApproach.coefs, rhs); }
  }
  void MatchEdges(Edge* a, Edge* b) {
    auto ab = Node(a->Center() - b->Center());
    auto b_right = make_unique<Cell>(*(a->left));
    auto a_right = make_unique<Cell>(*(b->left));
    a_right->Move(ab); b_right->Move(-ab);
    a->right = a_right.get(); b->right = b_right.get();
    ghost.emplace_back(move(a_right));
    ghost.emplace_back(move(b_right));
    auto dist = (a->left->Center() - a->right->Center()).norm();
    a->SetDist(dist); b->SetDist(dist);
  }
  // Data:
  vector<std::unique_ptr<Cell>> ghost;
  vector<Edge*> bdL, bdR, bdB, bdT;
  Real lower_[2], upper_[2];
};

template<int kOrder, class Physics>
struct OutFlow : public EdgeGroup<kOrder,Physics> {
  void PreProcess() override {
    for (auto& e : EdgeGroup<kOrder,Physics>::edge) {
      auto dist = (e->left->Center() - e->Center()).norm();
      e->SetDist(dist);
      e->right = e->left;
    }
  }
};

template<int kOrder, class Physics>
struct InFlow : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  using Solver = typename cfd::Solver<kOrder, Physics>;
  using Vr = VrApproach<kOrder, Physics>;
  using Mesh = typename Solver::MeshType;
  using State = Eigen::Matrix<Real, Physics::nEqual, 1>;
  using Coefs = Eigen::Matrix<Real, Mesh::nCoef * Physics::nEqual, 1>;
  // Construction:
  InFlow(std::function<void(Real, const Real*, Real*)>func) : Base(), func_(func) {}
  void CalculateBmats(Vr& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1]; WithoutDerivative(kOrder, distance, dp);
        auto c = e->left;
        e->Integrate([&](const Node& node) {
          return vr.GetMatAt(node.data(), *c, *c, normal, dp);
        }, &vr.B_mat[e->I()]);
    }
  }
  std::function<void(Real, const Real*, Real*)> func_;
  std::unordered_map<int, State> bdState;
  std::unordered_map<int, Coefs> bdCoef;
};

template<int kOrder, class Physics>
struct FarField : public EdgeGroup<kOrder,Physics> {
  // Types:
  using Base = EdgeGroup<kOrder,Physics>;
  using Vr = VrApproach<kOrder, Physics>;
  FarField(const Real* reference, int nVals) : Base() {
    refer_ = new Real[nVals];
    for (int i = 0; i < nVals; ++i) { refer_[i] = reference[i]; }
  }
  void CalculateBmats(Vr& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1]; WithoutDerivative(kOrder, distance, dp);
        auto c = e->left;
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
  using Vr = VrApproach<kOrder, Physics>;
  InviscWall() : Base() {}
  void CalculateBmats(Vr& vr) const override {
    for (auto& e : Base::edge) {
      Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
        Real dp[kOrder+1]; WithoutDerivative(kOrder, distance, dp);
        auto c = e->left;
        e->Integrate([&](const Node& node) {
          return vr.GetMatAt(node.data(), *c, *c, normal, dp);
        }, &vr.B_mat[e->I()]);
    }
  }
};
template<int kOrder, class Physics>
struct EdgeManager {
  // Type:
  using Solver = typename cfd::Solver<kOrder, Physics>;
  using Mesh = typename Solver::MeshType;
  using Group = EdgeGroup<kOrder,Physics>;
  using ConVar = typename Physics::ConVar;
  EdgeManager() = default;
  void PreProcess() { for (auto& [type, bd] : bdGroup) { bd->PreProcess(); } }
  void InitializeBndConds(DM dm, const BndConds* bc) {
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
        bdGroup[types[i]] = make_unique<InFlow<kOrder,Physics>>(bc->inflow);
        break;
      case BdCondType::OutFlow:
        bdGroup[types[i]] = make_unique<OutFlow<kOrder,Physics>>();
        break;
      case BdCondType::FarField:
        bdGroup[types[i]] = make_unique<FarField<kOrder,Physics>>(bc->refVal, 1);
        break;
      case BdCondType::InviscWall:
        bdGroup[types[i]] = make_unique<InviscWall<kOrder,Physics>>();
        break;
      default:
        break;
      }
    }
  }
  void ClassifyEdges(Mesh& mesh) {
    int               eStart, eEnd, type;
    DMPlexGetDepthStratum(mesh.dm, 1, &eStart, &eEnd); /* edges */
    for (int i = 0; i < mesh.NumLocalCells(); ++i) {
      for(int j = mesh.offset[i]; j < mesh.offset[i+1]; ++j) {
        int e = mesh.interface[j];
        if (mesh.edge[e]->right == nullptr) {
          DMGetLabelValue(mesh.dm, "Face Sets", e+eStart, &type);
          types.emplace(e, type);
          bdGroup[type]->edge.insert(mesh.edge[e].get());
        } else { // interior edges
          bdGroup[0]->edge.insert(mesh.edge[e].get());
        }
      }
    }
  }
  // Data:
  ConVar                                          cv;
  std::unordered_map<int, int>                    types;
  std::unordered_map<int, std::unique_ptr<Group>> bdGroup;
};

}  // cfd

#endif // INCLUDE_BNDCONDS_HPP_

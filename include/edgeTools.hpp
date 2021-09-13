/// @file edgeTools.hpp
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

#ifndef INCLUDE_EDGETOOLS_HPP_
#define INCLUDE_EDGETOOLS_HPP_

#include "defs.h"

namespace cfd {

template<int kOrder, class Physics>
struct Group {
  // Constants:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr int nEqual = Physics::nEqual;
  // Type:
  using Flux = typename Physics::Flux;
  using State = typename Physics::State;
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  // Virtual Functions:
  virtual void PreProcess() {
    for (auto& e : edge) {
      auto dist = (e->left->Center() - e->Center()).norm();
      e->SetDist(dist);
    }
  }
  virtual Matrix CalculateMat(const Edge<kOrder>* e, const Cell<kOrder>* left, 
                              const Cell<kOrder>* right) const {
    Matrix mat = Matrix::Zero();
    Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
    Real dp[kOrder+1]; InteriorDp(kOrder, distance, dp);
    e->Integrate([&](const Node& node) {
      return GetMatAt(node.data(), *left, *right, normal, dp);
    }, &mat);
    return mat;
  }
  virtual void UpdateRHS(const Array& priVar, const Array& coefs, Array& rhs) const {
    for (auto& e : edge) {
      auto flux = InteriorFlux(e, priVar, coefs, rhs);
      rhs.col(e->left->I())  -= flux;
      rhs.col(e->right->I()) += flux;
    }
  }
  // Static Functions:
  static constexpr Matrix GetMatAt(const Real* coord,
                                   const Cell<kOrder>& a, const Cell<kOrder>& b,
                                   const Real* normal, const Real* dp) {
    auto i = a.GetFuncTable(coord, normal);
    auto j = b.GetFuncTable(coord, normal);
    Matrix mat = Matrix::Zero();
    for (int m = 0; m != nCoef; ++m) {
      for (int n = 0; n != nCoef; ++n) {
        for (int k = 0; k != kOrder+1; ++k) mat(m,n) += dp[k] * i(n,k) * j(m,k);
      }
    }
    return mat;
  }
  static constexpr Flux InteriorFlux(const Edge<kOrder>* e, const Array& priVar,
                                     const Array& coefs, Array& rhs) {
    auto cell_l = e->left, cell_r = e->right;
    Flux flux_c = Flux::Zero();
    const Real normal[2] = {e->Nx(), e->Ny()};
    int i = cell_l->I(), j = cell_r->I();
    e->Integrate([&](const Node& p){
      State U_l = priVar.col(i) + coefs.block<nCoef, nEqual>(0, i*nEqual).transpose() *
                  cell_l->Functions(p.data());
      State U_r = priVar.col(j) + coefs.block<nCoef, nEqual>(0, j*nEqual).transpose() *
                  cell_r->Functions(p.data());
      return Physics::GetFlux(normal, U_l.data(), U_r.data());
    }, &flux_c);
    return flux_c;
  }
  EdgeSet<kOrder> edge;
};

template<int kOrder, class Physics>
struct Interior : public Group<kOrder,Physics> {
  // Types:
  using Base = Group<kOrder,Physics>;
  // Functions:
  Interior() : Group<kOrder,Physics>() {}
  void PreProcess() override {
    for (auto& e : Base::edge) {
      auto dist = (e->left->Center() - e->right->Center()).norm();
      e->SetDist(dist);
    }
  }
};

template<int kOrder, class Physics>
struct Periodic : public Group<kOrder,Physics> {
  // Types:
  using Base = Group<kOrder,Physics>;
  using Flux = typename Physics::Flux;
  using State = typename Physics::State;
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
    auto cmpY = [](Edge<kOrder>* a, Edge<kOrder>* b) {
      auto c_a = a->Center()(1), c_b = b->Center()(1);
      return c_a < c_b; };
    auto cmpX = [](Edge<kOrder>* a, Edge<kOrder>* b) {
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
  void UpdateRHS(const Array& priVar, const Array& coefs, Array& rhs) const override {
    for (auto& e : bdL) {
      auto flux = Base::InteriorFlux(e, priVar, coefs, rhs);
      rhs.col(e->left->I())  -= flux;
      rhs.col(e->right->I()) += flux;
    }
    for (auto& e : bdB) {
      auto flux = Base::InteriorFlux(e, priVar, coefs, rhs);
      rhs.col(e->left->I())  -= flux;
      rhs.col(e->right->I()) += flux;
    }
  }
  // Self Functions:
  void MatchEdges(Edge<kOrder>* a, Edge<kOrder>* b) {
    auto ab = Node(a->Center() - b->Center());
    auto b_right = make_unique<Cell<kOrder>>(*(a->left));
    auto a_right = make_unique<Cell<kOrder>>(*(b->left));
    a_right->Move(ab); b_right->Move(-ab);
    a->right = a_right.get(); b->right = b_right.get();
    ghost.emplace_back(move(a_right));
    ghost.emplace_back(move(b_right));
    auto dist = (a->left->Center() - a->right->Center()).norm();
    a->SetDist(dist); b->SetDist(dist);
  }
  // Data:
  vector<std::unique_ptr<Cell<kOrder>>> ghost;
  vector<Edge<kOrder>*> bdL, bdR, bdB, bdT;
  Real lower_[2], upper_[2];
};

template<int kOrder, class Physics>
struct OutFlow : public Group<kOrder,Physics> {
  // Type:
  using Flux = typename Physics::Flux;
  using State = typename Physics::State;
  // Constants:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr int nEqual = Physics::nEqual;
  // Functions:
  void UpdateRHS(const Array& priVar, const Array& coefs, Array& rhs) const override {
    for (auto& e : Group<kOrder,Physics>::edge) {
      auto cell = e->left; int i = cell->I();
      Flux flux_c = Flux::Zero();
      const Real normal[2] = {e->Nx(), e->Ny()};
      e->Integrate([&](const Node& p){
        State U = priVar.col(i) + coefs.block<nCoef, nEqual>(0, i*nEqual).transpose() *
                  cell->Functions(p.data());
        return Physics::GetFlux(normal, U.data());
      }, &flux_c);
      rhs.col(i) -= flux_c;
    }
  }
};

template<int kOrder, class Physics>
struct InFlow : public Group<kOrder,Physics> {
  // Constants:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  // Types:
  using Base = Group<kOrder,Physics>;
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  // Functions:
  InFlow(function<void(Real, const Real*, Real*)>func) : Base(), func_(func) {}
  Matrix CalculateMat(const Edge<kOrder>* e, const Cell<kOrder>* left, 
                      const Cell<kOrder>* right) const override {
    Matrix mat = Matrix::Zero();
    Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
    Real dp[kOrder+1]; WithoutDerivative(kOrder, distance, dp);
    e->Integrate([&](const Node& node) {
      return Base::GetMatAt(node.data(), *left, *right, normal, dp);
    }, &mat);
    return mat;
  }
  function<void(Real, const Real*, Real*)> func_;
};

template<int kOrder, class Physics>
struct FarField : public Group<kOrder,Physics> {
  // Constants:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  // Types:
  using Base = Group<kOrder,Physics>;
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  // Functions:
  FarField(const Real* reference, int nVals) : Base() {
    refer_ = new Real[nVals];
    for (int i = 0; i < nVals; ++i) { refer_[i] = reference[i]; }
  }
  ~FarField() { delete[] refer_;}
  virtual Matrix CalculateMat(const Edge<kOrder>* e, const Cell<kOrder>* left, 
                              const Cell<kOrder>* right) const {
    Matrix mat = Matrix::Zero();
    Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
    Real dp[kOrder+1]; WithoutDerivative(kOrder, distance, dp);
    e->Integrate([&](const Node& node) {
      return Base::GetMatAt(node.data(), *left, *right, normal, dp);
    }, &mat);
    return mat;
  }
  Real* refer_;
};

template<int kOrder, class Physics>
struct InviscWall : public Group<kOrder,Physics> {
  // Constants:
  static constexpr int nCoef = (kOrder+1)*(kOrder+2)/2-1; /**< Dofs -1 */
  static constexpr int nEqual = Physics::nEqual;
  // Types:
  using Base = Group<kOrder,Physics>;
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  using Flux = typename Physics::Flux;
  using State = typename Physics::State;
  // Functions:
  Matrix CalculateMat(const Edge<kOrder>* e, const Cell<kOrder>* left, 
                      const Cell<kOrder>* right) const override {
    Matrix mat = Matrix::Zero();
    Real normal[2] = {e->Nx(), e->Ny()}; Real distance = e->Distance();
    Real dp[kOrder+1]; WithoutDerivative(kOrder, distance, dp);
    e->Integrate([&](const Node& node) {
      return Base::GetMatAt(node.data(), *left, *right, normal, dp);
    }, &mat);
    return mat;
  }
  void UpdateRHS(const Array& priVar, const Array& coefs, Array& rhs) const override {
    for (auto& e : Group<kOrder,Physics>::edge) {
      auto cell = e->left; int i = cell->I();
      Flux flux_c = Flux::Zero();
      const Real normal[2] = {e->Nx(), e->Ny()};
      e->Integrate([&](const Node& p){
        State U = priVar.col(i) + coefs.block<nCoef, nEqual>(0, i*nEqual).transpose() *
                  cell->Functions(p.data());
        return Physics::GetWallFlux(normal, U.data());
      }, &flux_c);
      rhs.col(i) -= flux_c;
    }
  }
};

}  // cfd

#endif // INCLUDE_EDGETOOLS_HPP_

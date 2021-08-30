/// @file element.hpp
///
/// Definition of the class related to element and to metrics.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_GEOMETRY_ELEMENT_HPP_
#define INCLUDE_GEOMETRY_ELEMENT_HPP_

#include "defs.h"
#include "geometry/quadrature.hpp"

namespace cfd {

template <int kOrder>
class Edge
{
 public:
  static constexpr int nQuad = int(kOrder/2)+1;
  static constexpr LineQuad<nQuad> quad = LineQuad<nQuad>();
  static constexpr LineQuad<1> quad_1 = LineQuad<1>();
  static constexpr LineQuad<2> quad_2 = LineQuad<2>();
  static constexpr LineQuad<3> quad_3 = LineQuad<3>();
  static constexpr LineQuad<4> quad_4 = LineQuad<4>();
  Edge(int id, const Node& head, const Node& tail) : id_(id), head_(head), tail_(tail){
    measure_ = (head_ - tail_).norm();
  }
  int I() const { return id_; }
  const Node& Head() const { return head_; }
  const Node& Tail() const { return tail_; }
  Real Nx() const { return (tail_(1) - head_(1)) / measure_; }
  Real Ny() const { return (head_(0) - tail_(0)) / measure_; }
  Node Center() const { return (head_ + tail_) * 0.5; }
  Real Measure() const { return measure_; }
  Real Distance() const { return distance_; }
  void SetDist(Real dist) { distance_ = dist; }

  template <class Value, class Integrand>
  void Integrate(Integrand&& integrand, Value* value) const {
    Node points[nQuad];
    Quadrature(nQuad, points);
    for (int i = 0; i < nQuad; ++i) {
      *value += integrand(points[i]) * quad.w[i];
    }
    *value *= 0.5 * Measure();
  }
  void Quadrature(int npoints, Node* points) const {
    auto mid = Center();
    auto len = Tail() - Head();
    switch(npoints) {
      case 1:
        *points = mid;
        break;
      case 2:
        points[0] = mid + len * quad_2.x[0] * 0.5;
        points[1] = mid + len * quad_2.x[1] * 0.5;
        break;
      case 3:
        points[0] = mid + len * quad_3.x[0] * 0.5;
        points[1] = mid + len * quad_3.x[1] * 0.5;
        points[2] = mid + len * quad_3.x[2] * 0.5;
        break;
      case 4:
        points[0] = mid + len * quad_4.x[0] * 0.5;
        points[1] = mid + len * quad_4.x[1] * 0.5;
        points[2] = mid + len * quad_4.x[2] * 0.5;
        points[3] = mid + len * quad_4.x[3] * 0.5;
        break;
      default :
        break;
    }
  }
  /* Cells on both sides */
  Cell<kOrder>* left{nullptr};
  Cell<kOrder>* right{nullptr};

 private:
  int id_;
  Real measure_;
  Real distance_{0};
  const Node& head_;
  const Node& tail_;
};

template <>
class Cell<1>
{
 public:
  // Types:
  static constexpr int numCoef = 2;
  using Vector = Eigen::Matrix<Real, 2, 1>;
  using BasisF = Eigen::Matrix<Real, 2, 2>;
  // For Triangle
  Cell(int id, const Node& a, const Node& b, const Node& c) : id_(id) {
    measure_ = GetMeasure(a, b, c);
    center_ = (a + b + c) / 3;
    dx_inv_ = GetDelta(a(0), b(0), c(0));
    dy_inv_ = GetDelta(a(1), b(1), c(1));
  }
  // For Quadrangle
  Cell(int id, const Node& a, const Node& b, const Node& c, const Node& d) : id_(id) {
    measure_ = GetMeasure(a, b, c, d);
    center_ = (a + b + c + d) / 4;
    dx_inv_ = GetDelta(a(0), b(0), c(0), d(0));
    dy_inv_ = GetDelta(a(1), b(1), c(1), d(1));
  }
  int I() const {return id_; }
  Real Measure() const { return measure_; }
  void Move(const Node& begToEnd) { center_ += begToEnd; }
  // Accessors:
  static constexpr int Order() { return 1; }
  virtual int nCorner() const {}
  Real DxInv() const { return dx_inv_; }
  Real DyInv() const { return dy_inv_; }
  const Node& Center() const { return center_; }
  // For Triangle
  static Real GetMeasure(const Node& a, const Node& b, const Node& c) {
    auto measure = ((a(0) - b(0)) * (a(1) + b(1)) +
                    (b(0) - c(0)) * (b(1) + c(1)) +
                    (c(0) - a(0)) * (c(1) + a(1))) * 0.5;
    return Abs(measure);
  }
  static Real GetDelta(Real a, Real b, Real c) {
    auto delta = (max(max(a, b), c) - min(min(a, b), c)) * 0.5;
    return 1.0 / delta;
  }
  // For Quadrangle
  static Real GetMeasure(const Node& a, const Node& b, const Node& c, const Node& d) {
    auto measure = ((a(0) - c(0)) * (b(1) - d(1)) +
                    (d(0) - b(0)) * (a(1) - c(1))) * 0.5;
    return Abs(measure);
  }
  static Real GetDelta(Real a, Real b, Real c, Real d) {
    auto delta = (max(max(max(a, b), c), d) - min(min(min(a, b), c), d)) * 0.5;
    return 1.0 / delta;
  }
  // Basis Functions:
  Real F_0_0_0(const Real* coord) const { return (coord[0] - Center()(0)) * DxInv(); }
  Real F_0_1_0(const Real* coord) const { return DxInv(); }
  Real F_1_0_0(const Real* coord) const { return (coord[1] - Center()(1)) * DyInv(); }
  Real F_1_0_1(const Real* coord) const { return DyInv(); }
  // Normal Functions:
  Real F_0_N_1(const Real* coord, const Real* n) const { return F_0_1_0(coord) * n[0]; }
  Real F_1_N_1(const Real* coord, const Real* n) const { return F_1_0_1(coord) * n[1]; }
  Vector Functions(const Real* coord) const {
    Vector result;
    result << F_0_0_0(coord), F_1_0_0(coord);
    return result;
  }
  BasisF GetFuncTable(const Real* coord, const Real* normal) const {
    BasisF mat = BasisF::Zero();
    mat(0, 0) = F_0_0_0(coord), mat(0, 1) = F_0_N_1(coord, normal);
    mat(1, 0) = F_1_0_0(coord), mat(1, 1) = F_1_N_1(coord, normal);
    return mat;
  }

 private:
  int id_;
  Real measure_;
  Real dx_inv_, dy_inv_;
  Node center_;
};

template <>
class Cell<2> : public Cell<1>
{
 public:
  // Types:
  static constexpr int numCoef = 5;
  using Vector = Eigen::Matrix<Real, 5, 1>;
  using BasisF = Eigen::Matrix<Real, 5, 3>;
  Cell(int id, const Node& a, const Node& b, const Node& c) : Cell<1>(id, a, b, c) {
    xx_ = IntegrateTri([&](Real* coord) { return Pow(F_0_0_0(coord), 2);}, a, b, c) / Measure();
    xy_ = IntegrateTri([&](Real* coord) { return F_0_0_0(coord)*F_1_0_0(coord);}, a, b, c) / Measure();
    yy_ = IntegrateTri([&](Real* coord) { return Pow(F_1_0_0(coord), 2);}, a, b, c) / Measure();
  }
  // For Quadrangle
  Cell(int id, const Node& a, const Node& b, const Node& c, const Node& d) : Cell<1>(id, a, b, c, d) {
    xx_ = IntegrateQua([&](Real* coord) { return Pow(F_0_0_0(coord), 2);}, a, b, c, d) / Measure();
    xy_ = IntegrateQua([&](Real* coord) { return F_0_0_0(coord)*F_1_0_0(coord);}, a, b, c, d) / Measure();
    yy_ = IntegrateQua([&](Real* coord) { return Pow(F_1_0_0(coord), 2);}, a, b, c, d) / Measure();
  }
  static constexpr int Order() { return 2; }
  Real XX() const { return xx_; }
  Real XY() const { return xy_; }
  Real YY() const { return yy_; }

  // Basis Functions:
  Real F_2_0_0(const Real* coord) const { return Pow(F_0_0_0(coord), 2) - XX(); }
  Real F_2_1_0(const Real* coord) const { return F_0_0_0(coord) * DxInv() * 2; }
  Real F_2_2_0(const Real* coord) const { return Pow(DxInv(), 2) * 2; }
  // Normal Functions:
  Real F_2_N_1(const Real* coord, const Real* n) const { return F_2_1_0(coord) * n[0]; }
  Real F_2_N_2(const Real* coord, const Real* n) const { return F_2_2_0(coord) * n[0] * n[0]; }

  Real F_3_0_0(const Real* coord) const { return F_0_0_0(coord) * F_1_0_0(coord) - XY(); }
  Real F_3_1_0(const Real* coord) const { return DxInv() * F_1_0_0(coord); }
  Real F_3_0_1(const Real* coord) const { return F_0_0_0(coord) * DyInv(); }
  Real F_3_1_1(const Real* coord) const { return DxInv() * DyInv(); }
  // Normal Functions:
  Real F_3_N_1(const Real* coord, const Real* n) const { return F_3_1_0(coord) * n[0] + F_3_0_1(coord) * n[1]; }
  Real F_3_N_2(const Real* coord, const Real* n) const { return F_3_1_1(coord) * n[0] * n[1] * 2; }

  Real F_4_0_0(const Real* coord) const { return Pow(F_1_0_0(coord), 2) - YY(); }
  Real F_4_0_1(const Real* coord) const { return F_1_0_0(coord) * DyInv() * 2; }
  Real F_4_0_2(const Real* coord) const { return Pow(DyInv(), 2) * 2; }
  // Normal Functions:
  Real F_4_N_1(const Real* coord, const Real* n) const { return F_4_0_1(coord) * n[1]; }
  Real F_4_N_2(const Real* coord, const Real* n) const { return F_4_0_2(coord) * n[1] * n[1]; }

  Vector Functions(const Real* coord) const {
    Vector result;
    result(0) = F_0_0_0(coord); result(1) = F_1_0_0(coord);
    result(2) = Pow(result(0), 2) - XX();
    result(3) = result(0) * result(1) - XY();
    result(4) = Pow(result(1), 2) - YY();
    return result;
  }
  BasisF GetFuncTable(const Real* coord, const Real* normal) const {
    BasisF mat = BasisF::Zero();
    mat(0, 0) = F_0_0_0(coord), mat(0, 1) = F_0_N_1(coord, normal);
    mat(1, 0) = F_1_0_0(coord), mat(1, 1) = F_1_N_1(coord, normal);
    mat(2, 0) = F_2_0_0(coord), mat(2, 1) = F_2_N_1(coord, normal), mat(2, 2) = F_2_N_2(coord, normal);
    mat(3, 0) = F_3_0_0(coord), mat(3, 1) = F_3_N_1(coord, normal), mat(3, 2) = F_3_N_2(coord, normal);
    mat(4, 0) = F_4_0_0(coord), mat(4, 1) = F_4_N_1(coord, normal), mat(4, 2) = F_4_N_2(coord, normal);
    return mat;
  }
  template<class Integrand>
  static Real IntegrateTri(Integrand&& func,
      const Node& a, const Node& b, const Node& c)
  {
    Real integrand{0};
    Eigen::Matrix<Real, 3, 3> transform_mat;
    transform_mat << 1.0, 1.0, 1.0, a(0), b(0), c(0), a(1), b(1), c(1);
    for (int i = 0; i < 3; ++i) {
      Eigen::Matrix<Real, 3, 1> abc{tri_.a[i], tri_.b[i], tri_.c[i]};
      Eigen::Matrix<Real, 3, 1> xy = transform_mat * abc;
      integrand += func(xy.data()+1) * tri_.w[i];
    }
    return integrand *= GetMeasure(a, b, c);
  }

  template<class Integrand>
  static Real IntegrateQua(Integrand&& func,
      const Node& a, const Node& b, const Node& c, const Node& d)
  {
    Real integrand{0};
    for (int i = 0; i < 4; ++i) {
      auto N = GetN(qua_.a[i], qua_.b[i]);
      Real coord[2];
      coord[0] = N[0] * a(0) + N[1] * b(0)+ N[2] * c(0) + N[3] * d(0);
      coord[1] = N[0] * a(1) + N[1] * b(1)+ N[2] * c(1) + N[3] * d(1);
      integrand += func(coord) * qua_.w[i];
    }
    return integrand *= GetMeasure(a, b, c, d) * 0.25;
  }
  static Eigen::Matrix<Real, 4, 1> GetN(Real a, Real b) {
    auto N = Eigen::Matrix<Real, 4, 1>();
    N << (1 + a) * (1 + b), (1 - a) * (1 + b),
         (1 - a) * (1 - b), (1 + a) * (1 - b);
    return N * 0.25;
  }
 private:
  static constexpr TriQuad<2> tri_ = TriQuad<2>();
  static constexpr QuaQuad<2> qua_ = QuaQuad<2>();
  Real xx_, xy_, yy_;
};
template <>
class Cell<3> : public Cell<2>
{
 public:
  // Types:
  static constexpr int numCoef = 9;
  using Vector = Eigen::Matrix<Real, 9, 1>;
  using BasisF = Eigen::Matrix<Real, 9, 4>;
  Cell(int id, const Node& a, const Node& b, const Node& c) : Cell<2>(id, a, b, c) {
    xxx_ = IntegrateTri([&](const Real* coord) { return Pow(F_0_0_0(coord), 3);}, a, b, c) / Measure();
    xxy_ = IntegrateTri([&](const Real* coord) { return Pow(F_0_0_0(coord), 2)*F_1_0_0(coord);}, a, b, c) / Measure();
    xyy_ = IntegrateTri([&](const Real* coord) { return F_0_0_0(coord)*Pow(F_1_0_0(coord), 2);}, a, b, c) / Measure();
    yyy_ = IntegrateTri([&](const Real* coord) { return Pow(F_1_0_0(coord), 3);}, a, b, c) / Measure();
  }
  // For Quadrangle
  Cell(int id, const Node& a, const Node& b, const Node& c, const Node& d) : Cell<2>(id, a, b, c, d) {
    xxx_ = IntegrateQua([&](const Real* coord) { return Pow(F_0_0_0(coord), 3);}, a, b, c, d) / Measure();
    xxy_ = IntegrateQua([&](const Real* coord) { return Pow(F_0_0_0(coord), 2)*F_1_0_0(coord);}, a, b, c, d) / Measure();
    xyy_ = IntegrateQua([&](const Real* coord) { return F_0_0_0(coord)*Pow(F_1_0_0(coord), 2);}, a, b, c, d) / Measure();
    yyy_ = IntegrateQua([&](const Real* coord) { return Pow(F_1_0_0(coord), 3);}, a, b, c, d) / Measure();
  }
  static constexpr int Order() { return 3; }
  Real XXX() const { return xxx_; }
  Real XXY() const { return xxy_; }
  Real XYY() const { return xyy_; }
  Real YYY() const { return yyy_; }
  // Basis Functions:
  Real F_5_0_0(const Real* coord) const { return Pow(F_0_0_0(coord), 3) - XXX(); }
  Real F_5_1_0(const Real* coord) const { return Pow(F_0_0_0(coord), 2) * DxInv() * 3; }
  Real F_5_2_0(const Real* coord) const { return F_0_0_0(coord) * Pow(DxInv(), 2) * 6; }
  Real F_5_3_0(const Real* coord) const { return Pow(DxInv(), 3) * 6; }
  // Normal Functions:
  Real F_5_N_1(const Real* coord, const Real* n) const { return F_5_1_0(coord) * n[0]; }
  Real F_5_N_2(const Real* coord, const Real* n) const { return F_5_2_0(coord) * Pow(n[0], 2); }
  Real F_5_N_3(const Real* coord, const Real* n) const { return F_5_3_0(coord) * Pow(n[0], 3); }

  Real F_6_0_0(const Real* coord) const { return Pow(F_0_0_0(coord), 2) * F_1_0_0(coord) - XXY(); }
  Real F_6_1_0(const Real* coord) const { return F_0_0_0(coord) * F_1_0_0(coord) * DxInv() * 2; }
  Real F_6_0_1(const Real* coord) const { return Pow(F_0_0_0(coord), 2) * DyInv(); }
  Real F_6_2_0(const Real* coord) const { return F_1_0_0(coord) * Pow(DxInv(), 2) * 2; }
  Real F_6_1_1(const Real* coord) const { return F_0_0_0(coord) * DxInv() * DyInv() * 2; }
  Real F_6_2_1(const Real* coord) const { return Pow(DxInv(), 2) * DyInv() * 2; }
  // Normal Functions:
  Real F_6_N_1(const Real* coord, const Real* n) const { return F_6_1_0(coord) * n[0] + F_6_0_1(coord) * n[1]; }
  Real F_6_N_2(const Real* coord, const Real* n) const { return F_6_2_0(coord) * Pow(n[0], 2) + F_6_1_1(coord) * n[0] * n[1] * 2; }
  Real F_6_N_3(const Real* coord, const Real* n) const { return F_6_2_1(coord) * Pow(n[0], 2) * n[1] * 3; }

  Real F_7_0_0(const Real* coord) const { return F_0_0_0(coord) * Pow(F_1_0_0(coord), 2) - XYY(); }
  Real F_7_1_0(const Real* coord) const { return Pow(F_1_0_0(coord), 2) * DxInv(); }
  Real F_7_0_1(const Real* coord) const { return F_0_0_0(coord) * F_1_0_0(coord) * DyInv() * 2; }
  Real F_7_0_2(const Real* coord) const { return F_0_0_0(coord) * Pow(DyInv(), 2) * 2; }
  Real F_7_1_1(const Real* coord) const { return F_1_0_0(coord) * DxInv() * DyInv() * 2; }
  Real F_7_1_2(const Real* coord) const { return DxInv() * Pow(DyInv(), 2) * 2; }
  // Normal Functions:
  Real F_7_N_1(const Real* coord, const Real* n) const { return F_7_1_0(coord) * n[0] + F_7_0_1(coord) * n[1]; }
  Real F_7_N_2(const Real* coord, const Real* n) const { return F_7_1_1(coord) * n[0] * n[1] * 2 + F_7_0_2(coord) * Pow(n[1], 2); }
  Real F_7_N_3(const Real* coord, const Real* n) const { return F_7_1_2(coord) * n[0] * Pow(n[1], 2) * 3; }

  Real F_8_0_0(const Real* coord) const { return Pow(F_1_0_0(coord), 3) - YYY(); }
  Real F_8_0_1(const Real* coord) const { return Pow(F_1_0_0(coord), 2) * DyInv() * 3; }
  Real F_8_0_2(const Real* coord) const { return F_1_0_0(coord) * Pow(DyInv(), 2) * 6; }
  Real F_8_0_3(const Real* coord) const { return Pow(DyInv(), 3) * 6; }
  // Normal Functions:
  Real F_8_N_1(const Real* coord, const Real* n) const { return F_8_0_1(coord) * n[1]; }
  Real F_8_N_2(const Real* coord, const Real* n) const { return F_8_0_2(coord) * Pow(n[1], 2); }
  Real F_8_N_3(const Real* coord, const Real* n) const { return F_8_0_3(coord) * Pow(n[1], 3); }

  Vector Functions(const Real* coord) const {
    Vector result;
    result(0) = F_0_0_0(coord); result(1) = F_1_0_0(coord);
    result(2) = Pow(result(0), 2);
    result(3) = result(0) * result(1);
    result(4) = Pow(result(1), 2);
    result(5) = result(0) * result(2) - XXX();
    result(6) = result(1) * result(2) - XXY();
    result(7) = result(0) * result(4) - XYY();
    result(8) = result(1) * result(4) - YYY();
    result(2) -= XX();
    result(3) -= XY();
    result(4) -= YY();
    return result;
  }
  BasisF GetFuncTable(const Real* coord, const Real* normal) const {
    BasisF mat = BasisF::Zero();
    // One Degree:
    mat(0, 0) = F_0_0_0(coord), mat(0, 1) = F_0_N_1(coord, normal);
    mat(1, 0) = F_1_0_0(coord), mat(1, 1) = F_1_N_1(coord, normal);
    // Two Degree:
    mat(2, 0) = F_2_0_0(coord), mat(2, 1) = F_2_N_1(coord, normal), mat(2, 2) = F_2_N_2(coord, normal);
    mat(3, 0) = F_3_0_0(coord), mat(3, 1) = F_3_N_1(coord, normal), mat(3, 2) = F_3_N_2(coord, normal);
    mat(4, 0) = F_4_0_0(coord), mat(4, 1) = F_4_N_1(coord, normal), mat(4, 2) = F_4_N_2(coord, normal);
    // Three Degree:
    mat(5, 0) = F_5_0_0(coord), mat(5, 1) = F_5_N_1(coord, normal), mat(5, 2) = F_5_N_2(coord, normal), mat(5, 3) = F_5_N_3(coord, normal);
    mat(6, 0) = F_6_0_0(coord), mat(6, 1) = F_6_N_1(coord, normal), mat(6, 2) = F_6_N_2(coord, normal), mat(6, 3) = F_6_N_3(coord, normal);
    mat(7, 0) = F_7_0_0(coord), mat(7, 1) = F_7_N_1(coord, normal), mat(7, 2) = F_7_N_2(coord, normal), mat(7, 3) = F_7_N_3(coord, normal);
    mat(8, 0) = F_8_0_0(coord), mat(8, 1) = F_8_N_1(coord, normal), mat(8, 2) = F_8_N_2(coord, normal), mat(8, 3) = F_8_N_3(coord, normal);
    return mat;
  }
  template<class Integrand>
  static Real IntegrateTri(Integrand&& func,
      const Node& a, const Node& b, const Node& c)
  {
    Real integrand{0};
    Eigen::Matrix<Real, 3, 3> transform_mat;
    transform_mat << 1.0, 1.0, 1.0, a(0), b(0), c(0), a(1), b(1), c(1);
    for (int i = 0; i < 4; ++i) {
      Eigen::Matrix<Real, 3, 1> abc{tri_.a[i], tri_.b[i], tri_.c[i]};
      Eigen::Matrix<Real, 3, 1> xy = transform_mat * abc;
      integrand += func(xy.data()+1) * tri_.w[i];
    }
    return integrand *= GetMeasure(a, b, c);
  }
  template<class Integrand>
  static Real IntegrateQua(Integrand&& func,
      const Node& a, const Node& b, const Node& c, const Node& d)
  {
    Real integrand{0};
    for (int i = 0; i < 4; ++i) {
      auto N = GetN(qua_.a[i], qua_.b[i]);
      Real coord[2];
      coord[0] = N[0] * a(0) + N[1] * b(0)+ N[2] * c(0) + N[3] * d(0);
      coord[1] = N[0] * a(1) + N[1] * b(1)+ N[2] * c(1) + N[3] * d(1);
      integrand += func(coord) * qua_.w[i];
    }
    return integrand *= GetMeasure(a, b, c, d) * 0.25;
  }
 private:
  static constexpr TriQuad<3> tri_ = TriQuad<3>();
  static constexpr QuaQuad<3> qua_ = QuaQuad<3>();
  Real xxx_, xxy_, xyy_, yyy_;
};

template <int kOrder>
class Triangle : public Cell<kOrder>
{
 public:
  // Constructors:
  Triangle(int id, const Node& a, const Node& b, const Node& c) : Cell<kOrder>{id, a, b, c} {}
  int nCorner() const override { return 3; }
  Triangle(const Triangle<kOrder>&) = default;            // override default copy constructor
  Triangle<kOrder>& operator = (const Triangle<kOrder>&); // and assignment operator
};

template <int kOrder>
class Quadrangle : public Cell<kOrder>
{
 public:
  // Constructors:
  Quadrangle(int id, const Node& a, const Node& b, const Node& c, const Node& d) :
      Cell<kOrder>{id, a, b, c, d} {}
  int nCorner() const override { return 4; }
  Quadrangle(const Quadrangle<kOrder>&) = default;            // override default copy constructor
  Quadrangle<kOrder>& operator = (const Quadrangle<kOrder>&); // and assignment operator
};

}  // cfd

#endif // INCLUDE_GEOMETRY_ELEMENT_HPP_
/// @file cell.cpp
///
/// Solution of 2-D Euler Equations
/// on Unstructured Triangular Grids.
///
//  Features:
//  ~~~~~~~~~
//  # unstructured finite-volume scheme of Variational Reconstruction
//  # triangular elements only
//  # ideal gas model (other models possible)
//  # vanLeer/AUSM FVS scheme, Barth and Jespersen's limiter
//  # explicit multistage time-stepping scheme (Runge-Kutta type)
//
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#include "geometry/element.hpp"
#include "gtest/gtest.h"

class TriangleTest : public ::testing::Test {
 protected:
  using T1 = Triangle<1>;
  using T2 = Triangle<2>;
  using T3 = Triangle<3>;
  int id{0};
  Node a{0.0, 0.0}, b{1.0, 0.0}, c{0.0, 2.0};
  Real eps{1e-6};
};
TEST_F(TriangleTest, OneDegree) {
  T1 triangle = T1(id, a, b, c);
  EXPECT_EQ(triangle.I(), 0);
  EXPECT_EQ(triangle.Order(), 1);
  EXPECT_EQ(triangle.nCorner(), 3);
  EXPECT_EQ(triangle.Center()(0) * 3, 1.0);
  EXPECT_EQ(triangle.Center()(1) * 3, 2.0);
  EXPECT_EQ(triangle.Measure(), 1.0);
  EXPECT_EQ(triangle.DxInv(), 2.0);
  EXPECT_EQ(triangle.DyInv(), 1.0);
  // One Degree:
  Real coord[2] = {3.0, 2.0};
  EXPECT_EQ(triangle.F_0_0_0(coord), (coord[0] - triangle.Center()(0)) * triangle.DxInv());
  EXPECT_EQ(triangle.F_0_1_0(coord), 2);
  EXPECT_EQ(triangle.F_1_0_0(coord), (coord[1] - triangle.Center()(1)) * triangle.DyInv());
  EXPECT_EQ(triangle.F_1_0_1(coord), 1);
  // Move:
  Node move{1.0, 2.0};
  Real coord_move[2] = {4.0, 4.0};
  auto triangle_move = triangle;
  triangle_move.Move(move);
  EXPECT_EQ(triangle_move.Center()(0), triangle.Center()(0)+move(0));
  EXPECT_EQ(triangle_move.Center()(1), triangle.Center()(1)+move(1));
  EXPECT_NEAR(triangle_move.F_0_0_0(coord_move), triangle.F_0_0_0(coord), eps);
  EXPECT_NEAR(triangle_move.F_0_1_0(coord_move), triangle.F_0_1_0(coord), eps);
  EXPECT_NEAR(triangle_move.F_1_0_0(coord_move), triangle.F_1_0_0(coord), eps);
  EXPECT_NEAR(triangle_move.F_1_0_1(coord_move), triangle.F_1_0_1(coord), eps);
}
TEST_F(TriangleTest, TwoDegree) {
  T2 triangle = T2(id, a, b, c);
  EXPECT_EQ(triangle.Order(), 2);
  // Two Degree:
  auto xx = triangle.IntegrateTri([&](const Real* coord){
    return triangle.F_2_0_0(coord);}, a, b, c);
  EXPECT_NEAR(xx, 0, eps);
  EXPECT_NEAR(triangle.F_2_1_0(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_2_2_0(triangle.Center().data()), 8, eps);
  auto xy = triangle.IntegrateTri([&](const Real* coord){
    return triangle.F_3_0_0(coord);}, a, b, c);
  EXPECT_NEAR(xy, 0, eps);
  EXPECT_NEAR(triangle.F_3_1_0(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_3_0_1(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_3_1_1(triangle.Center().data()), 2, eps);
  auto yy = triangle.IntegrateTri([&](const Real* coord){
    return triangle.F_4_0_0(coord);}, a, b, c);
  EXPECT_NEAR(yy, 0, eps);
  EXPECT_NEAR(triangle.F_4_0_1(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_4_0_2(triangle.Center().data()), 2, eps);
  // Move:
  Real coord[2] = {3.0, 2.0};
  Node move{1.0, 2.0};
  Real coord_move[2] = {4.0, 4.0};
  auto triangle_move = triangle;
  triangle_move.Move(move);
  EXPECT_EQ(triangle_move.Center()(0), triangle.Center()(0)+move(0));
  EXPECT_EQ(triangle_move.Center()(1), triangle.Center()(1)+move(1));
  auto func_move = triangle_move.Functions(coord_move);
  EXPECT_NEAR(func_move(0), triangle.F_0_0_0(coord), eps);
  EXPECT_NEAR(func_move(1), triangle.F_1_0_0(coord), eps);
  EXPECT_NEAR(func_move(2), triangle.F_2_0_0(coord), eps);
  EXPECT_NEAR(func_move(3), triangle.F_3_0_0(coord), eps);
  EXPECT_NEAR(func_move(4), triangle.F_4_0_0(coord), eps);
}
TEST_F(TriangleTest, ThreeDegree) {
  T3 triangle = T3(id, a, b, c);
  EXPECT_EQ(triangle.Order(), 3);
  auto xxx = triangle.IntegrateTri([&](const Real* coord){
    return triangle.F_5_0_0(coord);}, a, b, c);
  EXPECT_NEAR(xxx, 0, eps);
  EXPECT_NEAR(triangle.F_5_1_0(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_5_2_0(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_5_3_0(triangle.Center().data()), 48, eps);
  auto xxy = triangle.IntegrateTri([&](const Real* coord){
    return triangle.F_6_0_0(coord);}, a, b, c);
  EXPECT_NEAR(xxy, 0, eps);
  EXPECT_NEAR(triangle.F_6_1_0(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_6_2_0(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_6_0_1(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_6_1_1(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_6_2_1(triangle.Center().data()), 8, eps);
  auto xyy = triangle.IntegrateTri([&](const Real* coord){
    return triangle.F_7_0_0(coord);}, a, b, c);
  EXPECT_NEAR(xyy, 0, eps);
  EXPECT_NEAR(triangle.F_7_1_0(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_7_0_1(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_7_0_2(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_7_1_1(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_7_1_2(triangle.Center().data()), 4, eps);
  auto yyy = triangle.IntegrateTri([&](const Real* coord){
    return triangle.F_8_0_0(coord);}, a, b, c);
  EXPECT_NEAR(yyy, 0, eps);
  EXPECT_NEAR(triangle.F_8_0_1(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_8_0_2(triangle.Center().data()), 0, eps);
  EXPECT_NEAR(triangle.F_8_0_3(triangle.Center().data()), 6, eps);
  // Move:
  Real coord[2] = {3.0, 2.0};
  Node move{1.0, 2.0};
  Real coord_move[2] = {4.0, 4.0};
  auto triangle_move = triangle;
  triangle_move.Move(move);
  EXPECT_EQ(triangle_move.Center()(0), triangle.Center()(0)+move(0));
  EXPECT_EQ(triangle_move.Center()(1), triangle.Center()(1)+move(1));
  auto func_move = triangle_move.Functions(coord_move);
  EXPECT_NEAR(func_move(0), triangle.F_0_0_0(coord), eps);
  EXPECT_NEAR(func_move(1), triangle.F_1_0_0(coord), eps);
  EXPECT_NEAR(func_move(2), triangle.F_2_0_0(coord), eps);
  EXPECT_NEAR(func_move(3), triangle.F_3_0_0(coord), eps);
  EXPECT_NEAR(func_move(4), triangle.F_4_0_0(coord), eps);
  EXPECT_NEAR(func_move(5), triangle.F_5_0_0(coord), eps);
  EXPECT_NEAR(func_move(6), triangle.F_6_0_0(coord), eps);
  EXPECT_NEAR(func_move(7), triangle.F_7_0_0(coord), eps);
  EXPECT_NEAR(func_move(8), triangle.F_8_0_0(coord), eps);
}
class QuadrangleTest : public ::testing::Test {
 protected:
  using Q1 = Quadrangle<1>;
  using Q2 = Quadrangle<2>;
  using Q3 = Quadrangle<3>;
  int id{0};
  Node a{-1.0, 0.0}, b{0.0, -1.0}, c{1.0, 0.0}, d{0.0, 1.0};
  Real eps{1e-6};
};
TEST_F(QuadrangleTest, OneDegree) {
  Q1 quadrangle = Q1(id, a, b, c, d);
  Real center[2] = {0.0, 0.0};
  EXPECT_EQ(quadrangle.I(), 0);
  EXPECT_EQ(quadrangle.Order(), 1);
  EXPECT_EQ(quadrangle.nCorner(), 4);
  EXPECT_EQ(quadrangle.Center()(0), center[0]);
  EXPECT_EQ(quadrangle.Center()(1), center[1]);
  EXPECT_EQ(quadrangle.Measure(), 2.0);
  EXPECT_EQ(quadrangle.DxInv(), 1.0);
  EXPECT_EQ(quadrangle.DyInv(), 1.0);
  // One Degree:
  Real coord[2] = {3.0, 2.0};
  EXPECT_EQ(quadrangle.F_0_0_0(coord), (coord[0] - center[0]) * quadrangle.DxInv());
  EXPECT_EQ(quadrangle.F_0_1_0(coord), 1);
  EXPECT_EQ(quadrangle.F_1_0_0(coord), (coord[1] - center[1]) * quadrangle.DyInv());
  EXPECT_EQ(quadrangle.F_1_0_1(coord), 1);
  Node move{1.0, 2.0};
  Real coord_move[2] = {4.0, 4.0};
  auto quadrangle_move = quadrangle;
  quadrangle_move.Move(move);
  EXPECT_EQ(quadrangle_move.Center()(0), quadrangle.Center()(0)+move(0));
  EXPECT_EQ(quadrangle_move.Center()(1), quadrangle.Center()(1)+move(1));
  EXPECT_NEAR(quadrangle_move.F_0_0_0(coord_move), quadrangle.F_0_0_0(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_0_1_0(coord_move), quadrangle.F_0_1_0(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_1_0_0(coord_move), quadrangle.F_1_0_0(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_1_0_1(coord_move), quadrangle.F_1_0_1(coord), eps);
}
TEST_F(QuadrangleTest, TwoDegree) {
  Q2 quadrangle = Q2(id, a, b, c, d);
  Real center[2] = {0.0, 0.0};
  EXPECT_EQ(quadrangle.Order(), 2);
  Real vol = quadrangle.IntegrateQua([&](const Real* coord){ return 1.0; }, a, b, c, d);
  EXPECT_EQ(vol, quadrangle.Measure());
  Real integrand = quadrangle.IntegrateQua([&](const Real* coord){
      return coord[0] * coord[1]; }, a, b, c, d);
  EXPECT_EQ(integrand, 0);
  // Two Degree:
  auto xx = quadrangle.IntegrateQua([&](const Real* coord){
    return quadrangle.F_2_0_0(coord);}, a, b, c, d);
  EXPECT_NEAR(xx, 0, eps);
  EXPECT_NEAR(quadrangle.F_2_1_0(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_2_2_0(quadrangle.Center().data()), 2, eps);
  auto xy = quadrangle.IntegrateQua([&](const Real* coord){
    return quadrangle.F_3_0_0(coord);}, a, b, c, d);
  EXPECT_NEAR(xy, 0, eps);
  EXPECT_NEAR(quadrangle.F_3_1_0(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_3_0_1(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_3_1_1(quadrangle.Center().data()), 1, eps);
  auto yy = quadrangle.IntegrateQua([&](const Real* coord){
    return quadrangle.F_4_0_0(coord);}, a, b, c, d);
  EXPECT_NEAR(yy, 0, eps);
  EXPECT_NEAR(quadrangle.F_4_0_1(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_4_0_2(quadrangle.Center().data()), 2, eps);
  // Move:
  Real coord[2] = {3.0, 2.0};
  Node move{1.0, 2.0};
  Real coord_move[2] = {4.0, 4.0};
  auto quadrangle_move = quadrangle;
  quadrangle_move.Move(move);
  EXPECT_EQ(quadrangle_move.Center()(0), quadrangle.Center()(0)+move(0));
  EXPECT_EQ(quadrangle_move.Center()(1), quadrangle.Center()(1)+move(1));
  EXPECT_NEAR(quadrangle_move.F_2_1_0(coord_move), quadrangle.F_2_1_0(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_2_2_0(coord_move), quadrangle.F_2_2_0(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_3_1_0(coord_move), quadrangle.F_3_1_0(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_3_0_1(coord_move), quadrangle.F_3_0_1(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_3_1_1(coord_move), quadrangle.F_3_1_1(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_4_0_1(coord_move), quadrangle.F_4_0_1(coord), eps);
  EXPECT_NEAR(quadrangle_move.F_4_0_2(coord_move), quadrangle.F_4_0_2(coord), eps);
}
TEST_F(QuadrangleTest, ThreeDegree) {
  Q3 quadrangle = Q3(id, a, b, c, d);
  EXPECT_EQ(quadrangle.Order(), 3);
  Real integrand = quadrangle.IntegrateQua([&](const Real* coord){
      return Pow(coord[0], 2) * coord[1]; }, a, b, c, d);
  EXPECT_EQ(integrand, 0);
  integrand = quadrangle.IntegrateQua([&](const Real* coord){
      return Pow(coord[1], 2) * coord[0]; }, a, b, c, d);
  EXPECT_EQ(integrand, 0);
  auto xxx = quadrangle.IntegrateQua([&](const Real* coord){
    return quadrangle.F_5_0_0(coord);}, a, b, c, d);
  EXPECT_NEAR(xxx, 0, eps);
  EXPECT_NEAR(quadrangle.F_5_1_0(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_5_2_0(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_5_3_0(quadrangle.Center().data()), 6, eps);
  auto xxy = quadrangle.IntegrateQua([&](const Real* coord){
    return quadrangle.F_6_0_0(coord);}, a, b, c, d);
  EXPECT_NEAR(xxy, 0, eps);
  EXPECT_NEAR(quadrangle.F_6_1_0(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_6_2_0(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_6_0_1(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_6_1_1(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_6_2_1(quadrangle.Center().data()), 2, eps);
  auto xyy = quadrangle.IntegrateQua([&](const Real* coord){
    return quadrangle.F_7_0_0(coord);}, a, b, c, d);
  EXPECT_NEAR(xyy, 0, eps);
  EXPECT_NEAR(quadrangle.F_7_1_0(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_7_0_1(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_7_0_2(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_7_1_1(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_7_1_2(quadrangle.Center().data()), 2, eps);
  auto yyy = quadrangle.IntegrateQua([&](const Real* coord){
    return quadrangle.F_8_0_0(coord);}, a, b, c, d);
  EXPECT_NEAR(yyy, 0, eps);
  EXPECT_NEAR(quadrangle.F_8_0_1(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_8_0_2(quadrangle.Center().data()), 0, eps);
  EXPECT_NEAR(quadrangle.F_8_0_3(quadrangle.Center().data()), 6, eps);
  // Move:
  Real coord[2] = {3.0, 2.0};
  Node move{1.0, 2.0};
  Real coord_move[2] = {4.0, 4.0};
  auto quadrangle_move = quadrangle;
  quadrangle_move.Move(move);
  EXPECT_EQ(quadrangle_move.Center()(0), quadrangle.Center()(0)+move(0));
  EXPECT_EQ(quadrangle_move.Center()(1), quadrangle.Center()(1)+move(1));
  auto func_move = quadrangle_move.Functions(coord_move);
  EXPECT_NEAR(func_move(0), quadrangle.F_0_0_0(coord), eps);
  EXPECT_NEAR(func_move(1), quadrangle.F_1_0_0(coord), eps);
  EXPECT_NEAR(func_move(2), quadrangle.F_2_0_0(coord), eps);
  EXPECT_NEAR(func_move(3), quadrangle.F_3_0_0(coord), eps);
  EXPECT_NEAR(func_move(4), quadrangle.F_4_0_0(coord), eps);
  EXPECT_NEAR(func_move(5), quadrangle.F_5_0_0(coord), eps);
  EXPECT_NEAR(func_move(6), quadrangle.F_6_0_0(coord), eps);
  EXPECT_NEAR(func_move(7), quadrangle.F_7_0_0(coord), eps);
  EXPECT_NEAR(func_move(8), quadrangle.F_8_0_0(coord), eps);
}
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
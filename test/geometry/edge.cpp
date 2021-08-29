/// @file edge.cpp
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
static char help[] = "Test for edge class in geometry.\n";

#include "geometry/element.hpp"
#include "gtest/gtest.h"

namespace cfd {

class EdgeTest : public ::testing::Test {
 protected:
  using Edge_T = Edge<3>;
  Node head{-1.0, 0.0}, tail{1.0, 0.0};
  Real eps = 1e-8;
  int id{0};
};
TEST_F(EdgeTest, ConstructorAndMethod) {
  auto e = Edge_T(id, head, tail);
  EXPECT_EQ(e.I(), 0);
  EXPECT_EQ(e.Head()(0), head(0));
  EXPECT_EQ(e.Tail()(0), tail(0));
  EXPECT_EQ(e.Head()(1), head(1));
  EXPECT_EQ(e.Tail()(1), tail(1));
  EXPECT_EQ(e.Measure(), 2);
  EXPECT_EQ(e.Center()(0) * 2, head(0) + tail(0));
  EXPECT_EQ(e.Center()(1) * 2, head(1) + tail(1));
}
TEST_F(EdgeTest, Quadature) {
  auto e = Edge_T(id, head, tail);
  Node* nodes = new Node[2]();
  e.Quadrature(2, nodes);
  EXPECT_EQ(nodes[0](0), -0.577350269189626);
  EXPECT_EQ(nodes[0](1), 0);
  EXPECT_EQ(nodes[1](0), +0.577350269189626);
  EXPECT_EQ(nodes[1](1), 0);
  delete[] nodes;
  auto func_1 = [&](const Node& point) { return point(0); };
  auto integrand = 0.0;
  e.Integrate(func_1, &integrand);
  EXPECT_EQ(integrand, 0);
  auto func_2 = [&](const Node& point) { return std::pow(point(0), 2); };
  integrand = 0.0;
  e.Integrate(func_2, &integrand);
  EXPECT_NEAR(integrand, 2.0/3, eps);
  auto func_3 = [&](const Node& point) { return std::pow(point(0), 3); };
  integrand = 0.0;
  e.Integrate(func_3, &integrand);
  EXPECT_EQ(integrand, 0);
}
TEST_F(EdgeTest, GaussPoint) {
  auto e = Edge_T(id, head, tail);
  Node* nodes = new Node[3]();
  e.Quadrature(3, nodes);
  Real integrand{0};
  for (int i = 0; i < 3; ++i) {
    integrand += Pow(nodes[i](0), 5) * e.quad_3.w[i];
  }
  EXPECT_NEAR(integrand, 0.0, eps);
  delete[] nodes;
  nodes = new Node[4]();
  e.Quadrature(4, nodes);
  integrand = 0;
  for (int i = 0; i < 4; ++i) {
    integrand += Pow(nodes[i](0), 7) * e.quad_4.w[i];
  }
  EXPECT_NEAR(integrand, 0.0, eps);
  delete[] nodes;
}

}  // cfd

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
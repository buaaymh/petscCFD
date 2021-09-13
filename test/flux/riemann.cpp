/// @file riemann.cpp
///
/// Solution of 2-D Euler Equations
/// on Unstructured Triangular/Quadrangle Grids.
///
//  Features:
//  ~~~~~~~~~
//  # unstructured finite-volume scheme of Variational Reconstruction
//  # ideal gas model (other models possible)
//  # Roe FDS scheme, Venkatakrishnan's limiter
//  # explicit multistage time-stepping scheme (Runge-Kutta type)
//
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================
static char help[] = "Test for Euler class for Roe riemann solver.\n";

#include "euler.hpp"
#include "gtest/gtest.h"

namespace cfd {

class EulerTest : public ::testing::Test {
 protected:
  using Flux  = typename Euler::Flux;
  using State = typename Euler::State;
  int dim = 2;
  Real gamma  = Gamma;
  Real normal[2] = {1.0, 0};
  Real vvel_l{1.5}, vvel_r{2.5};
  static void CompareFlux(int n, const Real* lhs, const Real* rhs) {
    for (int i = 0; i < n; ++i) { EXPECT_DOUBLE_EQ(lhs[i], rhs[i]); }
  }
};
TEST_F(EulerTest, PrimToCons) {
  Real rho{0.1}, u{+0.2}, v{-0.2}, p{0.3};
  Real  primitive[4] = {rho, u, v, p};
  Real  conserved[4] = {rho, rho*u, rho*v, p/(gamma-1)+0.5*rho*(u*u+v*v)};
  Euler::PrimToCons(primitive, primitive);
  CompareFlux(4, primitive, conserved);
}
TEST_F(EulerTest, ConsToPrim) {
  Real rho{0.1}, u{+0.2}, v{-0.2}, p{0.3};
  Flux primitive;
  primitive << rho, u, v, p;
  Flux conserved;
  conserved << rho, rho*u, rho*v, p/(gamma-1)+0.5*rho*(u*u+v*v);
  Euler::ConsToPrim(conserved.data(), conserved.data());
  CompareFlux(4, primitive.data(), conserved.data());
}
TEST_F(EulerTest, Sod) {
  Real  left[] = {1.000, 0.0, vvel_l, 1.0};
  Real right[] = {0.125, 0.0, vvel_r, 0.1};
  Real   mid[] = {0.426319, +0.927453, vvel_l, 0.303130};
  auto roe = Euler::GetFlux(normal, left, right);
  auto cmp = Euler::GetFlux(normal, mid);
  cout << "Sod problem" << endl;
  printf("mass : %.4f == %.4f\n", roe(0), cmp(0));
  printf("xmom : %.4f == %.4f\n", roe(1), cmp(1));
  printf("ymom : %.4f == %.4f\n", roe(2), cmp(2));
  printf("ener : %.4f == %.4f\n", roe(3), cmp(3));
}
TEST_F(EulerTest, ShockCollision) {
  Real  left[] = {5.99924, 19.5975, vvel_l, 460.894};
  Real right[] = {5.99242, 6.19633, vvel_r, 46.0950};
  Real   mid[] = {5.99924, 19.5975, vvel_l, 460.894};
  auto roe = Euler::GetFlux(normal, left, right);
  auto cmp = Euler::GetFlux(normal, mid);
  cout << "ShockCollision problem" << endl;
  printf("mass : %.4f == %.4f\n", roe(0), cmp(0));
  printf("xmom : %.4f == %.4f\n", roe(1), cmp(1));
  printf("ymom : %.4f == %.4f\n", roe(2), cmp(2));
  printf("ener : %.4f == %.4f\n", roe(3), cmp(3));
}
TEST_F(EulerTest, BlastFromLeft) {
  Real  left[] = {1.0, 0.0, vvel_l, 1e+3};
  Real right[] = {1.0, 0.0, vvel_r, 1e-2};
  Real   mid[] = {0.575062, 19.59745, vvel_l, 460.8938};
  auto roe = Euler::GetFlux(normal, left, right);
  auto cmp = Euler::GetFlux(normal, mid);
  cout << "BlastFromLeft problem" << endl;
  printf("mass : %.4f == %.4f\n", roe(0), cmp(0));
  printf("xmom : %.4f == %.4f\n", roe(1), cmp(1));
  printf("ymom : %.4f == %.4f\n", roe(2), cmp(2));
  printf("ener : %.4f == %.4f\n", roe(3), cmp(3));
}
TEST_F(EulerTest, BlastFromRight) {
  Real  left[] = {1.0, 0.0, vvel_l, 1e-2};
  Real right[] = {1.0, 0.0, vvel_r, 1e+2};
  Real   mid[] = {0.575113, -6.196328, vvel_r, 46.09504};
  auto roe = Euler::GetFlux(normal, left, right);
  auto cmp = Euler::GetFlux(normal, mid);
  cout << "BlastFromRight problem" << endl;
  printf("mass : %.4f == %.4f\n", roe(0), cmp(0));
  printf("xmom : %.4f == %.4f\n", roe(1), cmp(1));
  printf("ymom : %.4f == %.4f\n", roe(2), cmp(2));
  printf("ener : %.4f == %.4f\n", roe(3), cmp(3));
}
TEST_F(EulerTest, AlmostVaccumed) {
  Real  left[] = {1.0, -2.0, vvel_l, 0.4};
  Real right[] = {1.0, +2.0, vvel_r, 0.4};
  Real   mid[] = {0.21852, 0.0, vvel_l, 0.001894};
  auto roe = Euler::GetFlux(normal, left, right);
  auto cmp = Euler::GetFlux(normal, mid);
  cout << "AlmostVaccumed problem" << endl;
  printf("mass : %.4f == %.4f\n", roe(0), cmp(0));
  printf("xmom : %.4f == %.4f\n", roe(1), cmp(1));
  printf("ymom : %.4f == %.4f\n", roe(2), cmp(2));
  printf("ener : %.4f == %.4f\n", roe(3), cmp(3));
}
TEST_F(EulerTest, Vaccumed) {
  Real  left[] = {1.0, -4.0, vvel_l, 0.4};
  Real right[] = {1.0, +4.0, vvel_r, 0.4};
  Real   mid[] = {0.0, 0.0, vvel_l, 0.0};
  auto roe = Euler::GetFlux(normal, left, right);
  auto cmp = Euler::GetFlux(normal, mid);
  cout << "Vaccumed problem" << endl;
  printf("mass : %.4f == %.4f\n", roe(0), cmp(0));
  printf("xmom : %.4f == %.4f\n", roe(1), cmp(1));
  printf("ymom : %.4f == %.4f\n", roe(2), cmp(2));
  printf("ener : %.4f == %.4f\n", roe(3), cmp(3));
}

}  // cfd

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

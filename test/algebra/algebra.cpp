/// @file algebra.cpp
//
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#include <iostream>
#include "defs.hpp"
#include "gtest/gtest.h"

class AlgebraTest : public ::testing::Test {
 protected:
  static constexpr int kOrder = 2;
  static constexpr int nCoef = 5;
  static constexpr int nEqual = 2;
  using AGBR = Algebra<nEqual, kOrder>;
  using Basis = Eigen::Matrix<Real, nCoef, 1>;
  using State = Eigen::Matrix<Real, nEqual, 1>;
  using Matrix = Eigen::Matrix<Real, nCoef, nCoef>;
  using ValCoefs = Eigen::Matrix<Real, nCoef*nEqual, 1>;
  using CoefsMat = Eigen::Matrix<Real, nEqual, nCoef>;
  ValCoefs valCoefs{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  Basis basis = {1.0, 2.0, 3.0, 4.0, 5.0};
};
TEST_F(AlgebraTest, Gradient) {
  CoefsMat coefsMat;
  coefsMat << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
  State grad_1 = AGBR::Grad(basis.data(), valCoefs.data());
  State grad_2 = coefsMat * basis;
  EXPECT_EQ(grad_1(0), grad_2(0));
  EXPECT_EQ(grad_1(1), grad_2(1));
}
TEST_F(AlgebraTest, Add) {
  AGBR::Add(5, basis.data(), basis.data());
  EXPECT_EQ(basis(0), 2.0);
  EXPECT_EQ(basis(1), 4.0);
  EXPECT_EQ(basis(2), 6.0);
  EXPECT_EQ(basis(3), 8.0);
  EXPECT_EQ(basis(4), 10.0);
  AGBR::Sub(5, basis.data(), basis.data());
  EXPECT_EQ(basis(0), 0.0);
  EXPECT_EQ(basis(1), 0.0);
  EXPECT_EQ(basis(2), 0.0);
  EXPECT_EQ(basis(3), 0.0);
  EXPECT_EQ(basis(4), 0.0);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
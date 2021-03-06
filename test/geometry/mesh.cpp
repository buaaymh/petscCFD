/// @file mesh.cpp
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
static char help[] = "Test for mesh class in geometry.\n";

#include "geometry/mesh.hpp"
#include "solver.hpp"
#include "gtest/gtest.h"

namespace cfd {
  
using std::string;

class TriMeshTest : public ::testing::Test {
 protected:
  Mesh<2>   mesh = Mesh<2>();
  string    meshfile = "tri_test.msh";
  Real      eps{1e-8};
  string    dir{TEST_DATA_DIR};
};
TEST_F(TriMeshTest, ReadMesh) {
  /*
     3 -- [4] -- 2
     |  (1)   /  |
    [3]   [2]   [1]
     |  /   (0)  |
     0 -- [0] -- 1
  */
  mesh.ReadMeshFile(dir+meshfile);
  EXPECT_EQ(mesh.CountNodes(), 4);
  EXPECT_EQ(mesh.CountEdges(), 5);
  EXPECT_EQ(mesh.CountCells(), 2);
}
TEST_F(TriMeshTest, ReadNode) {
  /*
     3 -- [4] -- 2
     |  (1)   /  |
    [3]   [2]   [1]
     |  /   (0)  |
     0 -- [0] -- 1
  */
  mesh.ReadMeshFile(dir+meshfile);
  EXPECT_EQ(mesh.node[0](0), 0);
  EXPECT_EQ(mesh.node[0](1), 0);
  EXPECT_EQ(mesh.node[1](0), 1);
  EXPECT_EQ(mesh.node[1](1), 0);
  EXPECT_EQ(mesh.node[2](0), 1);
  EXPECT_EQ(mesh.node[2](1), 1);
  EXPECT_EQ(mesh.node[3](0), 0);
  EXPECT_EQ(mesh.node[3](1), 1);
}
TEST_F(TriMeshTest, ReadEdge) {
  /*
     3 -- [4] -- 2
     |  (1)   /  |
    [3]   [2]   [1]
     |  /   (0)  |
     0 -- [0] -- 1
  */
  mesh.ReadMeshFile(dir+meshfile);
  EXPECT_EQ(mesh.edge[0]->Center()(0), 0.5);
  EXPECT_EQ(mesh.edge[0]->Center()(1), 0.0);
  EXPECT_EQ(mesh.edge[0]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[0]->left, mesh.cell[0].get());
  EXPECT_EQ(mesh.edge[0]->right, nullptr);
  EXPECT_EQ(mesh.edge[1]->Center()(0), 1.0);
  EXPECT_EQ(mesh.edge[1]->Center()(1), 0.5);
  EXPECT_EQ(mesh.edge[1]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[1]->left, mesh.cell[0].get());
  EXPECT_EQ(mesh.edge[1]->right, nullptr);
  EXPECT_EQ(mesh.edge[2]->Center()(0), 0.5);
  EXPECT_EQ(mesh.edge[2]->Center()(1), 0.5);
  EXPECT_NEAR(mesh.edge[2]->Measure(), Sqrt(2.0), eps);
  EXPECT_EQ(mesh.edge[2]->left, mesh.cell[0].get());
  EXPECT_EQ(mesh.edge[2]->right,mesh.cell[1].get());
  EXPECT_EQ(mesh.edge[3]->Center()(0), 0.0);
  EXPECT_EQ(mesh.edge[3]->Center()(1), 0.5);
  EXPECT_EQ(mesh.edge[3]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[3]->left,  mesh.cell[1].get());
  EXPECT_EQ(mesh.edge[3]->right,nullptr);
  EXPECT_EQ(mesh.edge[4]->Center()(0), 0.5);
  EXPECT_EQ(mesh.edge[4]->Center()(1), 1.0);
  EXPECT_EQ(mesh.edge[4]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[4]->left,  mesh.cell[1].get());
  EXPECT_EQ(mesh.edge[4]->right,nullptr);
}
TEST_F(TriMeshTest, ReadCell) {
  /*
     3 -- [4] -- 2
     |  (1)   /  |
    [3]   [2]   [1]
     |  /   (0)  |
     0 -- [0] -- 1
  */
  mesh.ReadMeshFile(dir+meshfile);
  EXPECT_EQ(mesh.cell[0]->nCorner(), 3);
  EXPECT_EQ(mesh.cell[0]->Measure(), 0.5);
  EXPECT_NEAR(mesh.cell[0]->Center()(0), 2.0/3, eps);
  EXPECT_NEAR(mesh.cell[0]->Center()(1), 1.0/3, eps);
  EXPECT_EQ(mesh.cell[0]->Edge(0), 0);
  EXPECT_EQ(mesh.cell[0]->Edge(1), 1);
  EXPECT_EQ(mesh.cell[0]->Edge(2), 2);
  EXPECT_EQ(mesh.cell[1]->nCorner(), 3);
  EXPECT_EQ(mesh.cell[1]->Measure(), 0.5);
  EXPECT_NEAR(mesh.cell[1]->Center()(0), 1.0/3, eps);
  EXPECT_NEAR(mesh.cell[1]->Center()(1), 2.0/3, eps);
  EXPECT_EQ(mesh.cell[1]->Edge(0), 3);
  EXPECT_EQ(mesh.cell[1]->Edge(1), 2);
  EXPECT_EQ(mesh.cell[1]->Edge(2), 4);
}

class QuadMeshTest : public ::testing::Test {
 protected:
  Mesh<2> mesh = Mesh<2>();
  string      meshfile = "quad_test.msh";
  Real        eps{1e-8};
  string      dir{TEST_DATA_DIR};
};
TEST_F(QuadMeshTest, ReadMesh) {
  /*
     3 -- [2] -- 2
     |           |
    [3]   (0)   [1]
     |           |
     0 -- [0] -- 1
  */
  mesh.ReadMeshFile(dir+meshfile);
  EXPECT_EQ(mesh.CountNodes(), 4);
  EXPECT_EQ(mesh.CountEdges(), 4);
  EXPECT_EQ(mesh.CountCells(), 1);
}
TEST_F(QuadMeshTest, ReadEdge) {
  /*
     3 -- [2] -- 2
     |           |
    [3]   (0)   [1]
     |           |
     0 -- [0] -- 1
  */
  mesh.ReadMeshFile(dir+meshfile);
  EXPECT_EQ(mesh.edge[0]->Center()(0), 0.5);
  EXPECT_EQ(mesh.edge[0]->Center()(1), 0.0);
  EXPECT_EQ(mesh.edge[0]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[0]->left, mesh.cell[0].get());
  EXPECT_EQ(mesh.edge[0]->right, nullptr);
  EXPECT_EQ(mesh.edge[1]->Center()(0), 1.0);
  EXPECT_EQ(mesh.edge[1]->Center()(1), 0.5);
  EXPECT_EQ(mesh.edge[1]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[1]->left, mesh.cell[0].get());
  EXPECT_EQ(mesh.edge[1]->right, nullptr);
  EXPECT_EQ(mesh.edge[2]->Center()(0), 0.5);
  EXPECT_EQ(mesh.edge[2]->Center()(1), 1.0);
  EXPECT_EQ(mesh.edge[2]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[2]->left, mesh.cell[0].get());
  EXPECT_EQ(mesh.edge[2]->right,nullptr);
  EXPECT_EQ(mesh.edge[3]->Center()(0), 0.0);
  EXPECT_EQ(mesh.edge[3]->Center()(1), 0.5);
  EXPECT_EQ(mesh.edge[3]->Measure(), 1.0);
  EXPECT_EQ(mesh.edge[3]->left, mesh.cell[0].get());
  EXPECT_EQ(mesh.edge[3]->right,nullptr);
}
TEST_F(QuadMeshTest, ReadCell) {
  /*
     3 -- [2] -- 2
     |           |
    [3]   (0)   [1]
     |           |
     0 -- [0] -- 1
  */
  mesh.ReadMeshFile(dir+meshfile);
  EXPECT_EQ(mesh.cell[0]->nCorner(), 4);
  EXPECT_EQ(mesh.cell[0]->Measure(), 1);
  EXPECT_EQ(mesh.cell[0]->Center()(0), 0.5);
  EXPECT_EQ(mesh.cell[0]->Center()(1), 0.5);
  EXPECT_EQ(mesh.cell[0]->Edge(0), 0);
  EXPECT_EQ(mesh.cell[0]->Edge(1), 1);
  EXPECT_EQ(mesh.cell[0]->Edge(2), 2);
  EXPECT_EQ(mesh.cell[0]->Edge(3), 3);
}

}  // cfd

int main(int argc, char* argv[]) {
  PetscInitialize(&argc, &argv, (char*)0, help);
  ::testing::InitGoogleTest(&argc, argv);
  RUN_ALL_TESTS();
  PetscFinalize();
}

/// @file main.cpp
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
static char help[] = "2D Finite Volume Example.\n";

#include "solver.hpp"

namespace cfd {

/// Defination of BndConds to set boundary conditions in this demo
///
BndConds::BndConds() : lower{0.0, 0.0}, upper{1.0, 1.0} {}
BndConds::~BndConds() = default;

template<int kOrder>
class Advection {
 public:
  explicit Advection(const User* user) : user_(user) {};
  void Run()
  {
    const std::string dir{TEST_DATA_DIR};
    auto solver = Solver<kOrder, Linear>();
    solver.mesh.ReadMeshFile(dir + user_->filename);
    solver.SetupDataLayout();
    solver.SetBoundaryConditions(&bc_);
    solver.InitializeDS();
    solver.InitializeTS(user_);
    solver.InitializeSolution(User::InitFunc);
    solver.Calculate();
  }
 private:
  BndConds bc_;
  const User* user_;
};

}  // cfd

int main( int argc,char **args ) {

  using std::cout;
  using std::endl;
  
  cfd::User           user;
  PetscMPIInt         rank;

  // Initialize program
  PetscInitialize(&argc, &args, (char*)0, help);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0)
  {
    cout << endl
      << " *************************************************" << endl
      << " *                                               *" << endl
      << " *       2-D FLOW ON UNSTRUCTURED BOX MESH       *" << endl
      << " *                                               *" << endl
      << " *         RKFV FOR UNSTEADY EULER FLOWS         *" << endl
      << " *                                               *" << endl
      << " *     (c) Minghao Yang, CFD Solver project      *" << endl
      << " *                                               *" << endl
      << " *          Version 1.0 from 08/15/2021          *" << endl
      << " *                                               *" << endl
      << " *************************************************" << endl << endl;
  }
  if(user.order == 2) {
    auto model = cfd::Advection<2>(&user);
    model.Run();
  }
  PetscFinalize();
}

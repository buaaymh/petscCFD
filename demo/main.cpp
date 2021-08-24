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

using namespace std;

/// Main program of the flow solver.
///
/// @param argc   number of command line arguments
/// @param argv[] list of arguments
/// @return       EXIT_SUCCESS or an error code
///

template<int kOrder>
class Advection {
 public:
  explicit Advection(User* user) : user_(user) {};
  void Run()
  {
    const string dir{TEST_DATA_DIR};
    auto solver = Solver<kOrder, Linear>();
    solver.ReadMesh(dir+user_->meshfile);
    solver.SetupDataLayout();
    solver.ConstructMesh();
    solver.SetBoundaryConditions(user_);
    // solver.InitializeTS(user_);
    // solver.WriteSolution(user_);
  }
 protected:
  User* user_;
};

int main( int argc,char **args )
{
  User           user;
  PetscMPIInt    rank;
  int            order = 2;
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
  if(order == 2)
  {
    auto model = Advection<2>(&user);
    model.Run();
  }
  
  PetscFinalize();
}
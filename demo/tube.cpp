/// @file tube.cpp
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

#include "euler.hpp"
#include "solver.hpp"

// #define SGLPREC

namespace cfd {

struct User {
  int           order = 3;
  string        filename = "blast_200.msh";
  string        model = "blast";
  int           n_step = 500, output_interval = 2;
  Real          cfl = 1.0, tEnd = 0.038;

  static constexpr auto InitFunc = [](int dim, const Real* coord, int Nf, Real* u) {
    /***** Lax *****/
    // if (coord[0] < 0.5) {
    //   u[0] = 0.445; u[1] = 0.698; u[2] = 0.0; u[3] = 3.528;
    // } else {
    //   u[0] = 0.5;   u[1] = 0.0;   u[2] = 0.0; u[3] = 0.571;
    // }
    /***** Sod *****/
    // if (coord[0] < 0.5) {
    //   u[0] = 1.0;   u[1] = 0.0; u[2] = 0.0; u[3] = 1.0;
    // } else {
    //   u[0] = 0.125; u[1] = 0.0; u[2] = 0.0; u[3] = 0.1;
    // }
    /***** Blast *****/
    if (coord[0] < 0.1) {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0; u[3] = 1000.0;
    } else if (coord[0] > 0.9) {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0; u[3] = 100.0;
    } else {
      u[0] = 1.0; u[1] = 0.0; u[2] = 0.0; u[3] = 0.01;
    }
  };
};

/// Defination of BndConds to set boundary conditions in this demo
///
BndConds::BndConds() = default;
BndConds::~BndConds() = default;

template<int kOrder>
class Advection {
 public:
  explicit Advection(const User* user) : user_(user) {};
  void Run()
  {
    const std::string dir{TEST_DATA_DIR};
    auto solver = Solver<kOrder, Euler>();
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
  if(user.order == 1) {
    auto model = cfd::Advection<1>(&user);
    model.Run();
  } else if (user.order == 2) {
    auto model = cfd::Advection<2>(&user);
    model.Run();
  } else if (user.order == 3) {
    auto model = cfd::Advection<3>(&user);
    model.Run();
  }
  PetscFinalize();
}

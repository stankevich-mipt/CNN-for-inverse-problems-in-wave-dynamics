#pragma once

#include "Solvers/EulerSolver.h"
#include "Solvers/DormandPrinceSolver.h"
#include "Solvers/TrapezoidSolver.h"
#include "Solvers/RungeSolver4.h"
#include "Solvers/AdaptiveRKSolver.h"
#include "Solvers/RK2Solver.h"

// #include "Solvers/RungeSolver8_10.h"
// #include "Solvers/RungeSolver12_14.h"

#include "../NonlinearSolvers/Solvers/SimpleIterations.h"
#include "../NonlinearSolvers/Solvers/KrylovProjector.h"

template <typename Scalar, typename IndexType>
class SolversFactory
{
public:
  static DifferentialSolver<Scalar>* Build(const std::string& solverName)
  {
    if (solverName == "Euler")            return new EulerSolver<Scalar>();
    if (solverName == "DormandPrince")    return new DormandPrinceSolver<Scalar>();
    if (solverName == "Trapezoid")
    {
      NonlinearSystemSolver<Scalar, IndexType> *nonlinearSolver = new SimpleIterationsSolver<Scalar, IndexType>(0.1);
     // NonlinearSystemSolver<Scalar, IndexType> *nonlinearSolver = new KrylovProjectorSolver<Scalar, IndexType>(1.0, 100, 5);
      return new TrapezoidSolver<Scalar, IndexType>(nonlinearSolver);
    } 
    if (solverName == "RK4")              return new RungeSolver<Scalar>();
    if (solverName == "RKF45")            return new AdaptiveRungeSolver<Scalar>();
    if (solverName == "RK2")              return new RK2Solver<Scalar>();

    /*
    if (solverName == "RungeSolver8_10")  return new RungeSolver8_10<Scalar>();
    if (solverName == "RungeSolver12_14") return new RungeSolver12_14<Scalar>();
    */
    return nullptr;
  }
};

#pragma once

#include "../DifferentialSolvers/SolverState.h"
#include "NonlinearFunctionSystem.h"

template <typename Scalar, typename IndexType>
struct NonlinearSystemSolver
{
public:
  NonlinearSystemSolver(): func(nullptr), evalsCount(0)
  {
  }

  virtual ~NonlinearSystemSolver() = default;
  virtual void SetFunction(NonlinearFunctionSystem<Scalar, IndexType>* func) = 0;
  virtual bool AdvancePhase(Scalar* solution, Scalar& err) = 0;

  virtual void InitStep(const Scalar* initialGuess, Scalar* solution, const SolverState& solverState, IndexType maxIterationsCount, Scalar tolerance)
  {
    dimsCount = func->GetDimentionsCount(solverState);
    if (initialGuess != solution)
    {
      for(IndexType coordIndex = 0; coordIndex < dimsCount; ++coordIndex)
      {
        solution[coordIndex] = initialGuess[coordIndex];
      }
    }

    evalsCount = 0;
    err = tolerance * Scalar(2.0);

    this->solverState        = solverState;
    this->maxIterationsCount = maxIterationsCount;
    this->tolerance          = tolerance;
  }

  NonlinearFunctionSystem<Scalar, IndexType> *func;
  IndexType dimsCount;
  SolverState solverState;
  IndexType maxIterationsCount;
  Scalar tolerance;
  IndexType evalsCount;
  Scalar err;
};

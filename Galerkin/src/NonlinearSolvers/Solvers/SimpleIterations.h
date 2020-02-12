#pragma once

#include "../../DifferentialSolvers/SolverState.h"
#include "../NonlinearSystemSolver.h"

template <typename Scalar, typename IndexType>
struct SimpleIterationsSolver : public NonlinearSystemSolver<Scalar, IndexType>
{
  using NonlinearSystemSolver<Scalar, IndexType>::solverState;
  using NonlinearSystemSolver<Scalar, IndexType>::maxIterationsCount;
  using NonlinearSystemSolver<Scalar, IndexType>::tolerance;
  using NonlinearSystemSolver<Scalar, IndexType>::dimsCount;
  using NonlinearSystemSolver<Scalar, IndexType>::evalsCount;
  using NonlinearSystemSolver<Scalar, IndexType>::err;
  using NonlinearSystemSolver<Scalar, IndexType>::func;

  SimpleIterationsSolver(Scalar slopeMultiplier)
  {
    this->slopeMultiplier = slopeMultiplier;
  }

  virtual ~SimpleIterationsSolver()
  {
    delete [] nextEstimation;
  }

  virtual void SetFunction(NonlinearFunctionSystem<Scalar, IndexType> *func) override
  {
    this->func = func;
    dimsCount = func->GetMaxDimentionsCount();

    nextEstimation = new Scalar[dimsCount];
  }

  bool AdvancePhase(Scalar* solution, Scalar& err) override
  {
    ++evalsCount;
    err = func->GetValue(solution, solverState, nextEstimation);

    if (err < tolerance) return true;

    for(IndexType coordIndex = 0; coordIndex < dimsCount; ++coordIndex)
    {
      solution[coordIndex] += nextEstimation[coordIndex] * slopeMultiplier;
    }

    return false;
  }

  Scalar slopeMultiplier;

  Scalar *nextEstimation;
};

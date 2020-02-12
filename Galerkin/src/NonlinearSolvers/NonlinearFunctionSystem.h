#pragma once

#include "../DifferentialSolvers/SolverState.h"

template <typename Scalar, typename IndexType>
struct NonlinearFunctionSystem
{
  virtual ~NonlinearFunctionSystem() = default;

  virtual IndexType GetMaxDimentionsCount() const = 0;
  virtual IndexType GetDimentionsCount(const SolverState& solverState) const = 0;
  virtual Scalar GetValue(const Scalar* coords, const SolverState& solverState, Scalar* values) = 0;
};

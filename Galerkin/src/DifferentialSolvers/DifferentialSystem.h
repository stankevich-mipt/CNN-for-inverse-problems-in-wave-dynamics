#include "SolverState.h"

#pragma once
template <typename Scalar>
class DifferentialSystem
{
public:
  DifferentialSystem(int solverPhasesCount, int hierarchyLevelsCount = 1)
  {
    this->solverPhasesCount    = solverPhasesCount;
    this->hierarchyLevelsCount = hierarchyLevelsCount;
  }
  virtual ~DifferentialSystem() = default;

  virtual int GetDimentionsCount(const SolverState& solverState = SolverState()) const = 0;
  virtual int GetMaxDimentionsCount() const = 0;

  virtual void GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState = SolverState()) = 0;
  virtual Scalar GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, 
    const SolverState& solverState = SolverState(), const Scalar* mults = nullptr) = 0;

  virtual void GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState& solverState)
  {
  }

  virtual void GetCurrCoords(Scalar& time, Scalar* currCoords) const = 0;

  virtual void SetCurrCoords(Scalar time, const Scalar* newCoords, const Scalar* oldCoords, const SolverState& solverState)
  {
  }

  virtual void SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState)
  {
  }

  virtual void SetCurrCoords(Scalar time, const Scalar* oldCoords) = 0;

  int GetHierarchyLevelsCount() const
  {
    return hierarchyLevelsCount;
  }

  int GetMaxHierarchyLevel() const
  {
    return 1 << (hierarchyLevelsCount - 1);
  }

  int GetSolverPhasesCount() const
  {
    return solverPhasesCount;
  }

  virtual Scalar GetTimeStepPrediction() = 0;

private:
  int hierarchyLevelsCount;
  int solverPhasesCount;
};

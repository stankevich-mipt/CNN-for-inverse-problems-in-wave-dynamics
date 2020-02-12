#pragma once
#include "SolverState.h"
#include "DifferentialSystem.h"
#include <limits>

template<typename Scalar>
class DifferentialSolver
{
public:
  DifferentialSolver() = default;
  virtual ~DifferentialSolver() = default;

  virtual void    SetSystem(DifferentialSystem<Scalar>* system) = 0;
  virtual int     GetPhasesCount() const = 0;

  virtual void    InitStep(Scalar timeStep, Scalar tolerance, bool)
  {
    this->timeStep    = timeStep;
    this->tolerance   = tolerance;
    predictedStep = std::numeric_limits<Scalar>::max();
  }

  virtual void    InitStep(const SolverState&) = 0;
  virtual bool    AdvancePhase(const SolverState&) = 0;
  virtual void    AdvanceStep(const SolverState&) = 0;
  virtual void    RevertStep(Scalar currTime) = 0;

  virtual Scalar  GetLastStepError() const = 0;
  virtual Scalar  GetTimeStepPrediction() const = 0;

  Scalar GetCurrTime() const
  {
    return currTime;
  }
  Scalar GetCurrStep() const
  {
    return timeStep;
  }

protected:
  Scalar currTime;
  Scalar timeStep;
  Scalar tolerance;
  Scalar stepError;
  Scalar predictedStep;
};

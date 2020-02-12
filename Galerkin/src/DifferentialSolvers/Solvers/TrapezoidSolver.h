#include "../DifferentialSolver.h"
#include "../../NonlinearSolvers/NonlinearFunctionSystem.h"
#include "../../NonlinearSolvers/NonlinearSystemSolver.h"

template<typename Scalar, typename IndexType>
class TrapezoidSolver: public DifferentialSolver<Scalar>, public NonlinearFunctionSystem<Scalar, IndexType>
{
public:
  using DifferentialSolver<Scalar>::currTime;
  using DifferentialSolver<Scalar>::timeStep;
  using DifferentialSolver<Scalar>::tolerance;
  using DifferentialSolver<Scalar>::stepError;
  using DifferentialSolver<Scalar>::predictedStep;

  TrapezoidSolver(NonlinearSystemSolver<Scalar, IndexType> *nonlinearSolver)
  {
    this->nonlinearSolver = nonlinearSolver;
  }

  void SetSystem(DifferentialSystem<Scalar>* system) override
  {
    this->system = system;
    nonlinearSolver->SetFunction(this);

    currCoords  = new Scalar[system->GetMaxDimentionsCount()];
    nextCoords  = new Scalar[system->GetMaxDimentionsCount()];
    derivatives = new Scalar[system->GetMaxDimentionsCount()];
    probeCoords = new Scalar[system->GetMaxDimentionsCount()];

    tmp = new Scalar[system->GetMaxDimentionsCount()];

    oldCoords   = (system->GetHierarchyLevelsCount() > 1) ? new Scalar[system->GetMaxDimentionsCount()] : 0;

    predictedCoords = new Scalar[system->GetMaxDimentionsCount()];

    currTime = 0;
  }

  ~TrapezoidSolver()
  {
    delete [] currCoords;
    delete [] nextCoords;
    delete [] derivatives;
    delete [] probeCoords;

    if (system->GetHierarchyLevelsCount() > 1)
    {
      delete [] oldCoords;
    }

    delete [] predictedCoords;
  }

  int GetPhasesCount() const override
  {
    return 1;
  }

  void InitStep(Scalar timeStep, Scalar tolerance, bool updateInitialCoords) override
  {
    DifferentialSolver<Scalar>::InitStep(timeStep, tolerance, updateInitialCoords);
    if (system->GetHierarchyLevelsCount() == 1)
    {
      system->GetCurrCoords(currTime, currCoords);
      initialized = phaseAdvanced = false;
      iteration = 0;
    }
  }

  void InitStep(const SolverState& solverState) override
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      system->GetCurrCoords(currTime, currCoords, oldCoords, solverState);
      initialized = phaseAdvanced = false;
      iteration = 0;
    }
  }

  bool AdvancePhase(const SolverState& solverState) override
  {
    if (!initialized)
    {
      initialized = true;
      system->GetCurrDerivatives(derivatives, solverState);
      currStep = timeStep * (1 << solverState.hierarchyLevel);
      for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
      {
        nextCoords[coordIndex] = currCoords[coordIndex] + derivatives[coordIndex] * currStep;
      }
      nonlinearSolver->InitStep(nextCoords, nextCoords, solverState, 100, tolerance * currStep);
      system->SetCurrCoords(currTime + currStep, nextCoords, solverState);
    }

    if (!phaseAdvanced)
    {
      ++iteration;
      Scalar time;
      system->GetCurrCoords(time, nextCoords);

      /*
      for (int i = 0; i < system->GetDimentionsCount(solverState); ++i)
      {
        if (fabs(tmp[i] - nextCoords[i]) > 1e-8)
        {
          std::cout << "FUUU!!!\n";
        }
      }*/
      

      phaseAdvanced = nonlinearSolver->AdvancePhase(nextCoords, stepError);
      system->SetCurrCoords(currTime + currStep, nextCoords, solverState);

      stepError /= (tolerance * currStep);
      phaseAdvanced = phaseAdvanced && (stepError < Scalar(1.0) || iteration > nonlinearSolver->maxIterationsCount);
    }

    return phaseAdvanced;
  }

  void AdvanceStep(const SolverState& solverState) override
  {
    if (solverState.IsPreInitial())
    {
      currTime += timeStep;
    }
    if (system->GetHierarchyLevelsCount() > 1)
    {
      system->SetCurrCoords(currTime, nextCoords, oldCoords, solverState);
    } else
    {
      system->SetCurrCoords(currTime, nextCoords, solverState);
    }
  }

  void RevertStep(Scalar) override
  {
  }

  Scalar GetLastStepError() const override
  {
    return stepError;
  }

  Scalar GetTimeStepPrediction() const override
  {
    return timeStep;
  }

  IndexType GetMaxDimentionsCount() const override
  {
    return system->GetMaxDimentionsCount();
  }

  IndexType GetDimentionsCount(const SolverState& solverState) const override
  {
    return system->GetDimentionsCount(solverState);
  }

  virtual Scalar GetValue(const Scalar* coords, const SolverState& solverState, Scalar* values) override
  {
  //  system->SetCurrCoords(currTime + currStep, coords, solverState);
    system->GetCurrDerivatives(probeCoords, solverState);
    for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
    {
      predictedCoords[coordIndex] = currCoords[coordIndex] + currStep * Scalar(0.5) * (probeCoords[coordIndex] + derivatives[coordIndex]);
      values[coordIndex] = predictedCoords[coordIndex] - coords[coordIndex];
    }
    return system->GetErrorValue(currTime, coords, predictedCoords, solverState);
  }

private:
  Scalar* currCoords;
  Scalar* nextCoords;
  Scalar* probeCoords;
  Scalar* derivatives;
  Scalar* oldCoords;
  Scalar* predictedCoords;

  Scalar* tmp;

  Scalar currStep;

  bool phaseAdvanced;
  bool initialized;
  IndexType iteration;

  DifferentialSystem<Scalar>* system;
  NonlinearSystemSolver<Scalar, IndexType>* nonlinearSolver;
};

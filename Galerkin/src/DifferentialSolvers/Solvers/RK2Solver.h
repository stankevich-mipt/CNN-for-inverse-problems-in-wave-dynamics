#include "../DifferentialSolver.h"

// improved Euler method (Heun`s method)

template <typename Scalar>
class RK2Solver: public DifferentialSolver<Scalar>
{
public:
  using DifferentialSolver<Scalar>::timeStep;
  using DifferentialSolver<Scalar>::currTime;

  RK2Solver(): DifferentialSolver<Scalar>()
  {}

  void SetSystem(DifferentialSystem<Scalar>* system) override
  {
    this->system = system;
    currCoords = new Scalar[system->GetMaxDimentionsCount()];
    nextCoords = new Scalar[system->GetMaxDimentionsCount()];
    probeCoords = new Scalar[system->GetMaxDimentionsCount()];
    oldCoords = (system->GetHierarchyLevelsCount() > 1) ? new Scalar[system->GetMaxDimentionsCount()] : nullptr;

    k1 = new Scalar[system->GetMaxDimentionsCount()];
    k2 = new Scalar[system->GetMaxDimentionsCount()];

    currTime = Scalar(0.0);
  }
  ~RK2Solver()
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      delete[] oldCoords;
    }

    delete[] currCoords;
    delete[] nextCoords;
    delete[] probeCoords;

    delete[] k1;
    delete[] k2;
  }

  int GetPhasesCount() const override
  {
    return 2;
  }

  void InitStep(Scalar timeStep, Scalar tolerance, bool updateInitialCoords) override
  {
    DifferentialSolver<Scalar>::InitStep(timeStep, tolerance, updateInitialCoords);
    if (system->GetHierarchyLevelsCount() == 1)
    {
      system->GetCurrCoords(currTime, currCoords);
    }
  }

  void InitStep(const SolverState& solverState) override
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      system->GetCurrCoords(currTime, currCoords, oldCoords, solverState);
    }
  }

  bool AdvancePhase(const SolverState& solverState) override
  {
    Scalar currStep = timeStep * (1 << solverState.hierarchyLevel);

    switch (solverState.phaseIndex)
    {
      case 0:
      {
        system->GetCurrDerivatives(k1, solverState);

        #pragma omp parallel for
        for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + k1[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep, probeCoords, solverState);
      } break;

      case 1:
      {
        system->GetCurrDerivatives(k2, solverState);

        #pragma omp parallel for
        for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          nextCoords[coordIndex] = currCoords[coordIndex]
            + currStep * (k1[coordIndex] + k2[coordIndex]) * Scalar(0.5);
        }
      } break;
    }
    return true;
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
      system->SetCurrCoords(currTime, nextCoords);
    }
  }

  void RevertStep(Scalar) override
  {
    // never will be called
  }

  Scalar GetLastStepError() const override
  {
    return Scalar(1.0);
  }

  Scalar GetTimeStepPrediction() const override
  {
    return timeStep;
  }

private:
  Scalar *currCoords;
  Scalar *nextCoords;
  Scalar *probeCoords;
  Scalar *oldCoords;

  Scalar *k1;
  Scalar *k2;

  DifferentialSystem<Scalar> *system;
};

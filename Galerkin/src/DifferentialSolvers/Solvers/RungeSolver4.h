#include "../DifferentialSolver.h"

template <typename Scalar>
class RungeSolver: public DifferentialSolver<Scalar>
{
public:
  using DifferentialSolver<Scalar>::timeStep;
  using DifferentialSolver<Scalar>::currTime;

  RungeSolver(): DifferentialSolver<Scalar>()
  {}

  void SetSystem(DifferentialSystem<Scalar>* system) override
  {
    this->system = system;
    currCoords  = new Scalar[system->GetMaxDimentionsCount()];
    nextCoords  = new Scalar[system->GetMaxDimentionsCount()];
    probeCoords = new Scalar[system->GetMaxDimentionsCount()];
    oldCoords = (system->GetHierarchyLevelsCount() > 1) ? new Scalar[system->GetMaxDimentionsCount()] : nullptr;

    k1 = new Scalar[system->GetMaxDimentionsCount()];
    k2 = new Scalar[system->GetMaxDimentionsCount()];
    k3 = new Scalar[system->GetMaxDimentionsCount()];
    k4 = new Scalar[system->GetMaxDimentionsCount()];

    currTime = Scalar(0.0);
  }
  ~RungeSolver()
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      delete[] oldCoords;
    }

    delete [] currCoords;
    delete [] nextCoords;
    delete [] probeCoords;

    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
  }

  int GetPhasesCount() const override
  {
    return 4;
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
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(0.5) * k1[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(0.5), probeCoords, solverState);
      } break;

      case 1:
      {
        system->GetCurrDerivatives(k2, solverState);

        #pragma omp parallel for
        for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(0.5) * k2[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(0.5), probeCoords, solverState);
      } break;

      case 2:
      {
        system->GetCurrDerivatives(k3, solverState);

        #pragma omp parallel for
        for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + k3[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep, probeCoords, solverState);
      } break;

      case 3:
      {
        system->GetCurrDerivatives(k4, solverState);

        #pragma omp parallel for
        for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          nextCoords[coordIndex] = currCoords[coordIndex]
            + currStep * (k1[coordIndex] + Scalar(2.0) * k2[coordIndex] + Scalar(2.0) * k3[coordIndex] + k4[coordIndex]) / Scalar(6.0);
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
    return system->GetTimeStepPrediction() * Scalar(0.25);
  }

private:
  Scalar *currCoords;
  Scalar *nextCoords;
  Scalar *probeCoords;
  Scalar *oldCoords;

  Scalar *k1;
  Scalar *k2;
  Scalar *k3;
  Scalar *k4;

  DifferentialSystem<Scalar> *system;
};

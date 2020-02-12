#include "..//DifferentialSolver.h"

// Runge–Kutta–Fehlberg 45
template<typename Scalar>
class AdaptiveRungeSolver : public DifferentialSolver<Scalar>
{
public:
  using DifferentialSolver<Scalar>::timeStep;
  using DifferentialSolver<Scalar>::currTime;
  using DifferentialSolver<Scalar>::tolerance;
  using DifferentialSolver<Scalar>::stepError;
  using DifferentialSolver<Scalar>::predictedStep;

  virtual void SetSystem(DifferentialSystem<Scalar> *system) override
  {
    this->system = system;

    initialCoords = (system->GetHierarchyLevelsCount() > 1) ? new Scalar[system->GetMaxDimentionsCount()] : nullptr;
    oldCoords     = (system->GetHierarchyLevelsCount() > 1) ? new Scalar[system->GetMaxDimentionsCount()] : nullptr;

    currCoords  = new Scalar[system->GetMaxDimentionsCount()];
    nextCoords1 = new Scalar[system->GetMaxDimentionsCount()];
    nextCoords2 = new Scalar[system->GetMaxDimentionsCount()];
    probeCoords = new Scalar[system->GetMaxDimentionsCount()];

    k1 = new Scalar[system->GetMaxDimentionsCount()];
    k2 = new Scalar[system->GetMaxDimentionsCount()];
    k3 = new Scalar[system->GetMaxDimentionsCount()];
    k4 = new Scalar[system->GetMaxDimentionsCount()];
    k5 = new Scalar[system->GetMaxDimentionsCount()];
    k6 = new Scalar[system->GetMaxDimentionsCount()];
  }
  ~AdaptiveRungeSolver()
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      delete[] initialCoords;
      delete[] oldCoords;
    }

    delete [] currCoords;
    delete [] nextCoords1;
    delete [] nextCoords2;
    delete [] probeCoords;

    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
    delete [] k6;
  }

  virtual int GetPhasesCount() const override
  {
    return 6;
  }

  void InitStep(Scalar timeStep, Scalar tolerance, bool updateInitialCoords) override
  {
    DifferentialSolver<Scalar>::InitStep(timeStep, tolerance, updateInitialCoords);
    if (updateInitialCoords)
    {
      system->GetCurrCoords(currTime, system->GetHierarchyLevelsCount() > 1 ? initialCoords : currCoords);
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
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(0.25) * k1[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(0.25), probeCoords, solverState);
      } break;

      case 1:
      {
        system->GetCurrDerivatives(k2, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex]
            + Scalar(3.0 / 32.0) * k1[coordIndex] * currStep
            + Scalar(9.0 / 32.0) * k2[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(3.0 / 8.0), probeCoords, solverState);
      } break;

      case 2:
      {
        system->GetCurrDerivatives(k3, solverState);

        #pragma omp parallel for
        for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex]
            + Scalar(1932.0 / 2197.0) * k1[coordIndex] * currStep
            - Scalar(7200.0 / 2197.0) * k2[coordIndex] * currStep
            + Scalar(7296.0 / 2197.0) * k3[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(12.0 / 13.0), probeCoords, solverState);
      } break;

      case 3:
      {
        system->GetCurrDerivatives(k4, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex]
            + Scalar(439.0 / 216.0)  * k1[coordIndex] * currStep
            - Scalar(8.0)            * k2[coordIndex] * currStep
            + Scalar(3680.0 / 513.0) * k3[coordIndex] * currStep
            - Scalar(845.0 / 4104.0) * k4[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep, probeCoords, solverState);
      } break;

      case 4:
      {
        system->GetCurrDerivatives(k5, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex]
            - Scalar(8.0 / 27.0)      * k1[coordIndex] * currStep
            + Scalar(2.0)             * k2[coordIndex] * currStep
            - Scalar(3544.0 / 2565.0) * k3[coordIndex] * currStep
            + Scalar(1859.0 / 4104.0) * k4[coordIndex] * currStep
            - Scalar(11.0 / 40.0)     * k5[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(0.5), probeCoords, solverState);
      } break;
    
      case 5:
      {
        system->GetCurrDerivatives(k6, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          nextCoords1[coordIndex] = currCoords[coordIndex]
            + Scalar(25.0 / 216.0)     * k1[coordIndex] * currStep
            + Scalar(1408.0 / 2565.0)  * k3[coordIndex] * currStep
            + Scalar(2197.0 / 4104.0)  * k4[coordIndex] * currStep
            - Scalar(1.0 / 5.0)        * k5[coordIndex] * currStep;

          nextCoords2[coordIndex] = currCoords[coordIndex]
            + Scalar(16.0 / 135.0)     * k1[coordIndex] * currStep
            + Scalar(6656.0 / 12825.0) * k3[coordIndex] * currStep
            + Scalar(28561.0 / 56430.0)* k4[coordIndex] * currStep
            - Scalar(9.0 / 50.0)       * k5[coordIndex] * currStep
            + Scalar(2.0 / 55.0)       * k6[coordIndex] * currStep;
        } 

        stepError = (system->GetErrorValue(currTime, nextCoords1, nextCoords2, solverState) / tolerance) / currStep;
        predictedStep = timeStep * Scalar(pow(Scalar(1.0) / stepError, Scalar(1.0 / 4.0)));
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
      system->SetCurrCoords(currTime, nextCoords2, oldCoords, solverState);
    } else
    {
      system->SetCurrCoords(currTime, nextCoords2);
    }
  }

  void RevertStep(Scalar currTime) override
  {
    this->currTime = currTime;
    system->SetCurrCoords(currTime, system->GetHierarchyLevelsCount() > 1 ? initialCoords : currCoords);
  }

  Scalar GetLastStepError() const override
  {
    return stepError;
  }

  Scalar GetTimeStepPrediction() const override
  {
    return predictedStep;
  }

private:
  Scalar* initialCoords;
  Scalar* oldCoords;

  Scalar *currCoords;
  Scalar *nextCoords1;
  Scalar *nextCoords2;
  Scalar *probeCoords;

  Scalar *k1;
  Scalar *k2;
  Scalar *k3;
  Scalar *k4;
  Scalar *k5;
  Scalar *k6;

  DifferentialSystem<Scalar> *system;
};

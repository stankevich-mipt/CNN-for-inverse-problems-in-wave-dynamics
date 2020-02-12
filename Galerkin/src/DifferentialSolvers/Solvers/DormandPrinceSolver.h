#pragma once
#include "../DifferentialSolver.h"

template<typename Scalar>
class DormandPrinceSolver : public DifferentialSolver<Scalar>
{
public:
  using DifferentialSolver<Scalar>::timeStep;
  using DifferentialSolver<Scalar>::currTime;
  using DifferentialSolver<Scalar>::tolerance;
  using DifferentialSolver<Scalar>::stepError;
  using DifferentialSolver<Scalar>::predictedStep;

  DormandPrinceSolver() {}

  void SetSystem(DifferentialSystem<Scalar>* system) override
  {
    this->system  = system;
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
    k7 = new Scalar[system->GetMaxDimentionsCount()];

    currTime = Scalar(0.0);
  }

  virtual ~DormandPrinceSolver()
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      delete [] initialCoords;
      delete [] oldCoords;
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
    delete [] k7;
  }

  int GetPhasesCount() const override
  {
    return 7;
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
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(0.20) * k1[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(0.20), probeCoords, solverState);
      } break;

      case 1:
      {
        system->GetCurrDerivatives(k2, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(3.0 / 40.0) * k1[coordIndex] * currStep
                                                           + Scalar(9.0 / 40.0) * k2[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(3.0 / 10.0), probeCoords, solverState);
      } break;

      case 2:
      {
        system->GetCurrDerivatives(k3, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(44.0 / 45.0)  * k1[coordIndex] * currStep
                                                           - Scalar(56.0 / 15.0)  * k2[coordIndex] * currStep
                                                           + Scalar(32.0 / 9.0)   * k3[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(4.0 / 5.0), probeCoords, solverState);
      } break;

      case 3:
      {
        system->GetCurrDerivatives(k4, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(19372.0 / 6561.0) * k1[coordIndex] * currStep
                                                           - Scalar(25360.0 / 2187.0) * k2[coordIndex] * currStep
                                                           + Scalar(64448.0 / 6561.0) * k3[coordIndex] * currStep
                                                           - Scalar(212.0 / 729.0)    * k4[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep * Scalar(8.0 / 9.0), probeCoords, solverState);
      } break;

      case 4:
      {
        system->GetCurrDerivatives(k5, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(9017.0 / 3168.0)  * k1[coordIndex] * currStep
                                                           - Scalar(355.0 / 33.0)     * k2[coordIndex] * currStep
                                                           + Scalar(46732.0 / 5247.0) * k3[coordIndex] * currStep
                                                           + Scalar(49.0 / 176.0)     * k4[coordIndex] * currStep
                                                           - Scalar(5103.0 / 18656.0) * k5[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep, probeCoords, solverState);
      } break;

      case 5:
      {
        system->GetCurrDerivatives(k6, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          probeCoords[coordIndex] = currCoords[coordIndex] + Scalar(35.0 / 384.0)     * k1[coordIndex] * currStep
                                                           + Scalar(500.0 / 1113.0)   * k3[coordIndex] * currStep
                                                           + Scalar(125.0 / 192.0)    * k4[coordIndex] * currStep
                                                           - Scalar(2187.0 / 6784.0)  * k5[coordIndex] * currStep
                                                           + Scalar(11.0 / 84.0)      * k6[coordIndex] * currStep;
        }
        system->SetCurrCoords(currTime + currStep, probeCoords, solverState);
      } break;

      case 6:
      {
        system->GetCurrDerivatives(k7, solverState);

        #pragma omp parallel for
        for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
        {
          nextCoords1[coordIndex] = currCoords[coordIndex]
            + currStep * (  Scalar(35.0 / 384.0)        * k1[coordIndex]
                          + Scalar(500.0 / 1113.0)      * k3[coordIndex]
                          + Scalar(125.0 / 192.0)       * k4[coordIndex]
                          - Scalar(2187.0 / 6784.0)     * k5[coordIndex]
                          + Scalar(11.0 / 84.0)         * k6[coordIndex]);

          nextCoords2[coordIndex] = currCoords[coordIndex]
            + currStep * (  Scalar(5179.0 / 57600.0)    * k1[coordIndex]
                          + Scalar(7571.0 / 16695.0)    * k3[coordIndex]
                          + Scalar(393.0 / 640.0)       * k4[coordIndex]
                          - Scalar(92097.0 / 339200.0)  * k5[coordIndex]
                          + Scalar(187.0 / 2100.0)      * k6[coordIndex]
                          + Scalar(1.0 / 40.0)          * k7[coordIndex]);
        }
        stepError = (system->GetErrorValue(currTime, nextCoords1, nextCoords2, solverState) / tolerance) / currStep;
        predictedStep = timeStep * Scalar(pow(Scalar(1.0) / stepError, Scalar(1.0 / 5.0)));
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
  Scalar* currCoords;
  Scalar* oldCoords;
  Scalar* nextCoords1;
  Scalar* nextCoords2;
  Scalar* probeCoords;

  Scalar* k1;
  Scalar* k2;
  Scalar* k3;
  Scalar* k4;
  Scalar* k5;
  Scalar* k6;
  Scalar* k7;

  DifferentialSystem<Scalar>* system;
};

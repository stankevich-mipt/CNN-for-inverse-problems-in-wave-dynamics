#pragma once

#include <iostream>

struct SolverState
{
  SolverState(int hierarchyLevelsCount = 1, int phasesCount = 1): 
    globalStepIndex(0), hierarchyPhase(0), hierarchyLevel(0), phaseIndex(0)
  {
    this->hierarchyLevelsCount = hierarchyLevelsCount;
    this->phasesCount          = phasesCount;
  }

  void SetState(int globalStepIndex, int hierarchyPhase, int hierarchyLevel, int phaseIndex = 0)
  {
    this->globalStepIndex = globalStepIndex;
    this->hierarchyPhase  = hierarchyPhase;
    this->hierarchyLevel  = hierarchyLevel;
  }

  bool IsInitial() const
  {
    return IsFirstStep() && phaseIndex == 0;
  }

  bool IsFirstStep() const
  {
    return hierarchyPhase == 0 && hierarchyLevel == 0;
  }

  bool IsPreInitial() const
  {
    return IsLastStep() && phaseIndex + 1 == phasesCount;
  }

  bool IsLastStep() const
  {
    return hierarchyPhase == 1 && hierarchyLevel + 1 == hierarchyLevelsCount;
  }

  bool IsUseful() const
  {
    return hierarchyLevelsCount > 1 || hierarchyPhase == 1;
  }

  int Index() const
  {
    return GetLocalStepIndex() * hierarchyLevelsCount * 2 +
           hierarchyPhase      * hierarchyLevelsCount     +
           hierarchyLevel;
  }

  bool AllHierarchyLevelsCompleted() const
  {
    return IsInitial() && globalStepIndex % GetMaxHierarchyLevel() == 0; 
  }

  void SetInitialState(int globalStepIndex)
  {
    this->globalStepIndex = globalStepIndex;
    hierarchyPhase = hierarchyLevel = phaseIndex = 0;
  }

  SolverState GetNext(bool* allPhasesCompleted) const
  {
    SolverState nextState(*this);

    if (allPhasesCompleted)
    {
      *allPhasesCompleted = false;
    }

    nextState.phaseIndex++;
    if (nextState.phaseIndex == phasesCount)
    {
      if (allPhasesCompleted)
      {
        *allPhasesCompleted = true;
      }
      nextState.hierarchyLevel++;
      nextState.phaseIndex = 0;
    }

    if (nextState.hierarchyLevel == hierarchyLevelsCount)
    {
      nextState.hierarchyPhase++;
      nextState.hierarchyLevel = 0;
    }

    const int HierarchyPhasesCount = 2;
    if (nextState.hierarchyPhase == HierarchyPhasesCount)
    {
      nextState.globalStepIndex++;
      nextState.SetInitialState(nextState.globalStepIndex);
    }
    return nextState;
  }

  void SetPhaseIndex(int phaseIndex)
  {
    this->phaseIndex = phaseIndex;
  }

  void Print() const
  {
    std::cout << "globalStepIndex = " << globalStepIndex << "; " 
              << "hierarchyLevel = " << hierarchyLevel << "; " 
              << "hierarchyPhase = " << hierarchyPhase << "; " 
              << "phaseIndex = "     << phaseIndex << ";" << std::endl;
  }

  bool operator<(const SolverState& other) const
  {
    if (GetLocalStepIndex() < other.GetLocalStepIndex()) return GetLocalStepIndex() < other.GetLocalStepIndex();
    if (hierarchyPhase != other.hierarchyPhase) return hierarchyPhase < other.hierarchyPhase;
    if (hierarchyLevel != other.hierarchyLevel) return hierarchyLevel < other.hierarchyLevel;
    return phaseIndex < other.phaseIndex;
  }

  inline int GetMaxHierarchyLevel() const
  {
    return 1 << (hierarchyLevelsCount - 1);
  }

  int GetLocalStepIndex() const
  {
    return globalStepIndex % GetMaxHierarchyLevel();
  }

  int globalStepIndex;
  int hierarchyPhase;
  int hierarchyLevel;
  int phaseIndex;

  int hierarchyLevelsCount;
  int phasesCount;
};

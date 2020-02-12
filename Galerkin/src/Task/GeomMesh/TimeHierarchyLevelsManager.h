#pragma once

#include "../../DifferentialSolvers/SolverState.h"
#include "../../IO/Vtk/MeshVtkWriter.h"
#include "BucketStorage.h"
#include "GeomMesh/GeomMesh.h"

#include <queue>
#include <algorithm>
#include <sstream>
#include <assert.h>
#include <omp.h>

template <typename Space>
class TimeHierarchyLevelsManager
{
public:
  SPACE_TYPEDEFS

  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<int> IntTypeCellInfo;

  void Initialize(IndexType cellsCount, 
    IndexType hierarchyLevelsCount, IndexType solverPhasesCount,
    bool externalInitialization);

  IndexType GetLevel(IndexType cellIndex) const;
  void      SetLevel(IndexType cellIndex, IndexType hierarchyLevel);

  bool NeedToUpdate(IndexType cellIndex, const SolverState& solverState, bool* auxCell = nullptr) const;
  bool UseHalfStepSolution(IndexType cellIndex, const SolverState&, bool auxCell, bool toRead) const;
  bool UseHalfStepSolutionForNeighbour(IndexType cellIndex, const SolverState&, bool auxCell, IndexType correspondingCellIndex) const;

  void BuildTimeHierarchyLevels(const GeomMesh<Space>* const geomMesh, 
    const std::vector<Scalar>& cellMaxWaveSpeeds, bool allowCollisions, Scalar minTimeStep);
  void SaveToVtk(const GeomMesh<Space>* const geomMesh, IndexType globalStepIndex) const;

  IndexType GetComputationalTotalCost() const;
  IndexType GetComputationalTotalCost(const std::vector<IndexType>& cellsDomainIds, IndexType domainsCount) const;
  IndexType GetComputationalCost(IndexType cellIndex) const;

private:
  void BuildNeigbourCellsInfo(const GeomMesh<Space>* const geomMesh);

  BucketStorage<IndexType> timeHierarchyLevels;
  std::vector<char>        isCellUsed;
  std::vector<IndexType>   distToCell;
  IndexType                maxNeighbourCellsCount;
  std::vector<IndexType>   neighbourCells;
  IndexType                cellsCount;
  IndexType                hierarchyLevelsCount;
  IndexType                solverPhasesCount;
  IndexType                threadsCount;

  struct CellNeighboursInfo
  {
    // it`s about time hierarchy
    CellNeighboursInfo(): equal(false), higher(false), lower(false), basicEqual(false)
    {
    }
    bool equal;
    bool higher;
    bool lower;
    bool basicEqual;
  };

  std::vector<CellNeighboursInfo> cellNeighboursInfos;
};

template <typename Space>
void TimeHierarchyLevelsManager<Space>::Initialize(IndexType cellsCount, IndexType hierarchyLevelsCount, 
  IndexType solverPhasesCount, bool externalInitialization)
{
  this->cellsCount           = cellsCount;
  this->hierarchyLevelsCount = hierarchyLevelsCount;
  this->solverPhasesCount    = solverPhasesCount;

  maxNeighbourCellsCount = 1;

  IndexType perimeterCells = 1;
  for (IndexType phaseIndex = 0; phaseIndex < solverPhasesCount; ++phaseIndex)
  {
    perimeterCells *= (phaseIndex == 0) ? Space::FacesPerCell : (Space::FacesPerCell - 1);
    maxNeighbourCellsCount += perimeterCells;
  }

  threadsCount = 1;
  #pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
    {
      threadsCount = omp_get_num_threads();
    }
  }

  cellNeighboursInfos.resize(cellsCount);
  isCellUsed.resize(threadsCount * cellsCount, 0);
  distToCell.resize(threadsCount * cellsCount, 0);
  neighbourCells.resize(threadsCount * maxNeighbourCellsCount);

  if (!externalInitialization)
  {
    timeHierarchyLevels.Initialize(hierarchyLevelsCount, cellsCount);
  }
}

template <typename Space>
inline typename Space::IndexType TimeHierarchyLevelsManager<Space>::GetLevel(IndexType cellIndex) const
{
  return timeHierarchyLevels[cellIndex];
}

template <typename Space>
void TimeHierarchyLevelsManager<Space>::SetLevel(IndexType cellIndex, IndexType hierarchyLevel)
{
  timeHierarchyLevels.UpdateBucket(cellIndex, hierarchyLevel);
}

template <typename Space>
bool TimeHierarchyLevelsManager<Space>::NeedToUpdate(IndexType cellIndex, const SolverState& solverState, bool* auxCell) const
{
  bool needToUpdate = false;

  switch (solverState.hierarchyPhase)
  {
    case 0:
    {
      if ((cellNeighboursInfos[cellIndex].lower || cellNeighboursInfos[cellIndex].basicEqual) && 
        timeHierarchyLevels[cellIndex] == solverState.hierarchyLevel + 1 &&
        solverState.globalStepIndex % (1 << (solverState.hierarchyLevel + 1)) == 0)
      {
        if (auxCell) *auxCell = false;
        needToUpdate = true;
      }

      if (cellNeighboursInfos[cellIndex].higher && 
        timeHierarchyLevels[cellIndex] == solverState.hierarchyLevel + 0 &&
        solverState.globalStepIndex % (1 << (solverState.hierarchyLevel + 1)) == 0)
      {
        if (auxCell) *auxCell = true;
        needToUpdate = true;
      }
    } break;
    case 1:
    {
      bool updateCurrentCell = timeHierarchyLevels[cellIndex] == solverState.hierarchyLevel &&
        solverState.globalStepIndex % (1 << solverState.hierarchyLevel) == 0;

      bool updateCorrespondingCell = false;

      updateCorrespondingCell |= 
        cellNeighboursInfos[cellIndex].equal && 
        timeHierarchyLevels[cellIndex] + 0 == solverState.hierarchyLevel && 
        solverState.globalStepIndex % (1 << solverState.hierarchyLevel) == 0;

      updateCorrespondingCell |= 
        cellNeighboursInfos[cellIndex].lower && 
        timeHierarchyLevels[cellIndex] - 1 == solverState.hierarchyLevel && 
        solverState.globalStepIndex % (1 << solverState.hierarchyLevel) == 0;

      updateCorrespondingCell |= 
        cellNeighboursInfos[cellIndex].higher && 
        timeHierarchyLevels[cellIndex] + 1 == solverState.hierarchyLevel && 
        solverState.globalStepIndex % (1 << solverState.hierarchyLevel) == 0;

      if (updateCurrentCell || updateCorrespondingCell)
      {
        if (auxCell) *auxCell = !updateCurrentCell;
        needToUpdate = true;
      }
    } break;
  }

  return needToUpdate;
}

template <typename Space>
inline bool TimeHierarchyLevelsManager<Space>::UseHalfStepSolution(
  IndexType cellIndex, const SolverState& solverState, bool auxCell, bool toRead) const
{
  // it`s known that we need to update this cell
  if (GetLevel(cellIndex) == 0) return false;
  bool isInitial = solverState.phaseIndex == 0 && solverState.globalStepIndex % (1 << GetLevel(cellIndex)) == 0;

  bool updateCurrentCell = GetLevel(cellIndex) == solverState.hierarchyLevel && 
    solverState.globalStepIndex % (1 << solverState.hierarchyLevel) == 0;

  return !updateCurrentCell && (!isInitial || !toRead) && (auxCell == (solverState.hierarchyPhase == 1));
}

template <typename Space>
bool TimeHierarchyLevelsManager<Space>::UseHalfStepSolutionForNeighbour(
  IndexType cellIndex, const SolverState& solverState, bool auxCell, IndexType correspondingCellIndex) const
{
  if (GetLevel(correspondingCellIndex) == 0) return false;

  bool useHalfStepSolutionForCurrentCell = UseHalfStepSolution(cellIndex, solverState, auxCell, true);

  bool auxCorrespondingCell;
  if (NeedToUpdate(correspondingCellIndex, solverState, &auxCorrespondingCell))
  {
    return UseHalfStepSolution(correspondingCellIndex, solverState, auxCorrespondingCell, true);
  }

  if (solverState.globalStepIndex % (1 << GetLevel(correspondingCellIndex)) == 0) return false;
  if (solverState.globalStepIndex % (1 << GetLevel(correspondingCellIndex)) ==
    (1 << (GetLevel(correspondingCellIndex) - 1))) return true;

  assert(0);
  return false;
}

template <typename Space>
void TimeHierarchyLevelsManager<Space>::BuildTimeHierarchyLevels(const GeomMesh<Space>* const geomMesh,
  const std::vector<Scalar>& cellMaxWaveSpeeds,  bool allowCollisions, Scalar minTimeStep)
{
  // compute time hierarchy levels
  for (IndexType cellIndex = 0; cellIndex < cellsCount; ++cellIndex)
  { 
    if (allowCollisions && geomMesh->IsBoundaryCell(cellIndex))
    {
      // border cells have 1st hierarchy level
      timeHierarchyLevels.UpdateBucket(cellIndex, 0);
      continue;
    }
    Scalar cellTimeStep = geomMesh->GetMinHeight(cellIndex) / cellMaxWaveSpeeds[cellIndex];
    IndexType currentLevel = timeHierarchyLevels[cellIndex];
    // TODO: check it
    IndexType newLevel = std::min(currentLevel, IndexType(Log2(cellTimeStep / minTimeStep)) );
    //IndexType newLevel = std::min(currentLevel, std::max<IndexType>(IndexType(cellTimeStep / minTimeStep), 1) - 1);

    timeHierarchyLevels.UpdateBucket(cellIndex, newLevel);
  }

  IndexType hierarchyLevel = 0;
  int step = 0; // TODO: unused variable?
  while (hierarchyLevel < hierarchyLevelsCount)
  {
    typedef typename BucketStorage<IndexType>::ElementIterator ElementIterator;

    if (hierarchyLevel > 0)
    {
      // boundary cells
      for (ElementIterator it = timeHierarchyLevels.Begin(hierarchyLevel - 1); it != timeHierarchyLevels.End(hierarchyLevel - 1); ++it)
      {
        IndexType cellIndex = *it;
        if (distToCell[cellIndex] == solverPhasesCount && timeHierarchyLevels[cellIndex] == hierarchyLevel - 1 + step)
        {
          for(IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
          {
            IndexType correspondingCellIndex = geomMesh->GetCorrespondingCellIndex(cellIndex, faceNumber);
            if (correspondingCellIndex != IndexType(-1) && timeHierarchyLevels[correspondingCellIndex] > IndexType(hierarchyLevel))
            {
              timeHierarchyLevels.UpdateBucket(correspondingCellIndex, hierarchyLevel);
            }
          }
        }
      }
    }

    std::queue<IndexType> cellsQueue;
    for (ElementIterator it = timeHierarchyLevels.Begin(hierarchyLevel); it != timeHierarchyLevels.End(hierarchyLevel); ++it)
    {
      IndexType cellIndex = *it;
      cellsQueue.push(cellIndex);
      isCellUsed[cellIndex] = 1;
    }

    std::fill(distToCell.begin(), distToCell.end(), 0);

    while (!cellsQueue.empty())
    {
      IndexType cellIndex = cellsQueue.front();
      cellsQueue.pop();

      for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
      {
        IndexType correspondingCellIndex = geomMesh->GetCorrespondingCellIndex(cellIndex, faceNumber);
        if (correspondingCellIndex != IndexType(-1) && 
            !isCellUsed[correspondingCellIndex] && 
            distToCell[cellIndex] < solverPhasesCount)
        {
          cellsQueue.push(correspondingCellIndex);
          isCellUsed[correspondingCellIndex] = 1;
          distToCell[correspondingCellIndex] = distToCell[cellIndex] + 1;
          timeHierarchyLevels.UpdateBucket(correspondingCellIndex, hierarchyLevel);
        }
      }
    }

    if (hierarchyLevel == 0) 
    {
      ++hierarchyLevel;
    } else
    {
      step++;
      if (step == 2)
      {
        ++hierarchyLevel;
        step = 0;
      }
    }
  }

  timeHierarchyLevels.Cache();
  BuildNeigbourCellsInfo(geomMesh);
}

template <typename Space>
void TimeHierarchyLevelsManager<Space>::SaveToVtk(const GeomMesh<Space>* const geomMesh, IndexType globalStepIndex) const
{
  std::vector<int> hierarchyLevels(cellsCount);
  for (IndexType cellIndex = 0; cellIndex < cellsCount; ++cellIndex)
  {
    hierarchyLevels[cellIndex] = int(timeHierarchyLevels[cellIndex]);
  }

  MeshVtkWriter< Space, IntTypeCellInfo > meshWriter;

  std::ostringstream fileName;
  fileName << "out/TimeHierarchyLevels[" << globalStepIndex << "].vtk";
  meshWriter.Write(fileName.str(), geomMesh->nodes, geomMesh->cells, hierarchyLevels);
}

template <typename Space>
typename TimeHierarchyLevelsManager<Space>::IndexType 
  TimeHierarchyLevelsManager<Space>::GetComputationalTotalCost(const std::vector<IndexType>& cellsDomainIds, IndexType domainsCount) const
{
  std::vector<IndexType> domainsCosts(domainsCount, 0);

  SolverState solverState(hierarchyLevelsCount, solverPhasesCount);
  IndexType maxHierarchyLevel = 1 << (hierarchyLevelsCount - 1);

  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(cellsCount); ++cellIndex)
  {
    IndexType cellCost = 0;
    for (IndexType globalStepIndex = 0; globalStepIndex < maxHierarchyLevel; ++globalStepIndex)
    {
      for (IndexType hierarchyPhase = 0; hierarchyPhase < 2; ++hierarchyPhase)
      {
        for (IndexType hierarchyLevel = 0; hierarchyLevel < hierarchyLevelsCount; ++hierarchyLevel)
        {
          solverState.SetState(globalStepIndex, hierarchyPhase, hierarchyLevel);
          if (NeedToUpdate(cellIndex, solverState))
          {
            ++cellCost;
          }
        }
      }
    }
    #pragma omp atomic
    domainsCosts[cellsDomainIds[cellIndex]] += cellCost;
  }

  IndexType domainMaxCost = 0;
  for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
  {
    domainMaxCost = std::max(domainMaxCost, domainsCosts[domainIndex]);
  }

  return domainMaxCost / maxHierarchyLevel;
}

template <typename Space>
typename TimeHierarchyLevelsManager<Space>::IndexType TimeHierarchyLevelsManager<Space>::GetComputationalTotalCost() const
{
  return GetComputationalTotalCost(std::vector<IndexType>(cellsCount), 1);
}

template <typename Space>
typename TimeHierarchyLevelsManager<Space>::IndexType TimeHierarchyLevelsManager<Space>::GetComputationalCost(IndexType cellIndex) const
{
  int maxHierarchyLevel = 1 << (int(hierarchyLevelsCount) - 1);

  int cellCost = 0;
  for (int globalStepIndex = 0; globalStepIndex < maxHierarchyLevel; ++globalStepIndex)
  {
    for (IndexType hierarchyPhase = 0; hierarchyPhase < 2; ++hierarchyPhase)
    {
      for (IndexType hierarchyLevel = 0; hierarchyLevel < hierarchyLevelsCount; ++hierarchyLevel)
      {
        SolverState solverState(hierarchyLevelsCount, solverPhasesCount);
        solverState.SetState(globalStepIndex, hierarchyPhase, hierarchyLevel);
        if (NeedToUpdate(cellIndex, solverState))
        {
          cellCost++;
        }
      }
    }
  }
  return cellCost;
}

template <typename Space>
void TimeHierarchyLevelsManager<Space>::BuildNeigbourCellsInfo(const GeomMesh<Space>* const geomMesh)
{
  std::fill(isCellUsed.begin(), isCellUsed.end(), 0);
  std::fill(distToCell.begin(), distToCell.end(), 0);

  // neighbour cells info
  for (int phase = 0; phase < 2; ++phase)
  {
    #pragma omp parallel
    {
      int threadIndex = omp_get_thread_num();
      int segmentBegin = (threadIndex + 0) * cellsCount / threadsCount;
      int segmentEnd   = (threadIndex + 1) * cellsCount / threadsCount;

      for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
      {
        if (phase == 0)
        {
          cellNeighboursInfos[cellIndex] = CellNeighboursInfo();
        }
        std::queue<IndexType> cellsQueue;
        int neighbourCellsCount = 0;

        cellsQueue.push(cellIndex);
        neighbourCells[threadIndex * maxNeighbourCellsCount + neighbourCellsCount++] = cellIndex;
        isCellUsed[threadIndex * cellsCount + cellIndex] = 1;

        while (!cellsQueue.empty())
        {
          IndexType neighbourCellIndex = cellsQueue.front();
          cellsQueue.pop();

          for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
          {
            IndexType correspondingCellIndex = geomMesh->GetCorrespondingCellIndex(neighbourCellIndex, faceNumber);
            if (correspondingCellIndex != IndexType(-1) && 
                !isCellUsed[threadIndex * cellsCount + correspondingCellIndex] && 
                distToCell[ threadIndex * cellsCount + neighbourCellIndex] < solverPhasesCount)
            {
              switch (phase)
              {
                case 0:
                {
                  cellNeighboursInfos[cellIndex].lower  |= timeHierarchyLevels[cellIndex] >  timeHierarchyLevels[correspondingCellIndex];
                  cellNeighboursInfos[cellIndex].equal  |= timeHierarchyLevels[cellIndex] == timeHierarchyLevels[correspondingCellIndex];
                  cellNeighboursInfos[cellIndex].higher |= timeHierarchyLevels[cellIndex] <  timeHierarchyLevels[correspondingCellIndex];
                } break;
                case 1:
                {
                  cellNeighboursInfos[correspondingCellIndex].basicEqual |= 
                    cellNeighboursInfos[cellIndex].lower && timeHierarchyLevels[cellIndex] == timeHierarchyLevels[correspondingCellIndex] &&
                    !cellNeighboursInfos[correspondingCellIndex].lower;

                  assert(!(cellNeighboursInfos[cellIndex].lower && timeHierarchyLevels[cellIndex] < timeHierarchyLevels[correspondingCellIndex]));
                } break;
              }
              cellsQueue.push(correspondingCellIndex);
              neighbourCells[threadIndex * maxNeighbourCellsCount + neighbourCellsCount++] = correspondingCellIndex;
              isCellUsed[threadIndex * cellsCount + correspondingCellIndex] = 1;
              distToCell[threadIndex * cellsCount + correspondingCellIndex] = distToCell[threadIndex * cellsCount + neighbourCellIndex] + 1;
            }
          }
        }
      
        for (int cellNumber = 0; cellNumber < neighbourCellsCount; ++cellNumber)
        {
          IndexType cellToUpdateIndex = neighbourCells[threadIndex * maxNeighbourCellsCount + cellNumber];
          isCellUsed[threadIndex * cellsCount + cellToUpdateIndex] = 0;
          distToCell[threadIndex * cellsCount + cellToUpdateIndex] = 0;
        }
      }
    }
  }
}

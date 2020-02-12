#include "../../VolumeMethod/FunctionGetters/ContactFinder.h"

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellSolution(IndexType cellIndex, Vector refCoords, bool halfStepSolution) const
{ 
  typename System::ValueType result(Scalar(0.0));

  for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refCoords, functionIndex);
    for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        (halfStepSolution ? halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
          cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
      result.values[valueIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellSolution(IndexType cellIndex, IndexType pointIndex, PointType pointType, bool halfStepSolution) const
{
  typename System::ValueType result(Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = 0;

    switch (pointType)
    {
      case Basis: 
        basisFunctionValue = basisPointFunctionValues[functionIndex][pointIndex]; 
       break;
      case CellNode:
        basisFunctionValue = cellNodeBasisFunctionValues[functionIndex][pointIndex];
      break;
      default:
        assert(0);
    }

    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        (halfStepSolution ? halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
          cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
      result.values[valueIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetCellSolution(IndexType cellIndex, Vector globalPoint, bool halfStepSolution) const
{
  Vector points[Space::NodesPerCell];
  GetCellVertices(cellIndex, points);
  Vector refCoords = GlobalToRefVolumeCoords(globalPoint, points);
  return GetRefCellSolution(cellIndex, refCoords, halfStepSolution);
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellSolution(Scalar* coeffs, Vector refCoords) const
{
  typename System::ValueType result;
  result.SetZeroValues();

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refCoords, functionIndex);
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient = coeffs[valueIndex * functionsCount + functionIndex];
      result.values[valueIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::MediumParameters VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellParams(Scalar* coeffs, Vector refCoords) const
{
  typename System::MediumParameters result;
  std::fill(result.params, result.params + MediumParameters::ParamsCount, Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refCoords, functionIndex);
    for (IndexType paramIndex = 0; paramIndex < MediumParameters::ParamsCount; ++paramIndex)
    {
      Scalar basisFunctionCoefficient = coeffs[paramIndex * functionsCount + functionIndex];
      result.params[paramIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetCellSolution(IndexType cellIndex, Scalar* coeffs, Vector globalPoint) const
{
  Vector points[Space::NodesPerCell];
  GetCellVertices(cellIndex, points);
  Vector refCoords = GlobalToRefVolumeCoords(globalPoint, points);
  return GetRefCellSolution(coeffs, refCoords);
}

template <typename Space, typename FunctionSpace, typename System>
int VolumeMeshCommon<Space, FunctionSpace, System>::
  GetDimentionsCount(const SolverState& solverState) const
{
  return hierarchyDimentionsCount[solverState.Index()];
}

template <typename Space, typename FunctionSpace, typename System>
int VolumeMeshCommon<Space, FunctionSpace, System>::
  GetMaxDimentionsCount() const
{
  return dimsCount * functionsCount * cells.size();
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::GetCurrCoords(Scalar& time, Scalar* currCoords) const
{
  time = this->time;
  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
  {
    for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
      {
        currCoords[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] = 
          cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];
      }
    }
  } 
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState& solverState)
{
  time = this->time;
  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];
    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    bool auxCell;
    bool useHalfStepSolution;

    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            currCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] = 
              (useHalfStepSolution ?
               halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]:
               cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
          }
        }

        useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            oldCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] = 
              (useHalfStepSolution ?
               halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]:
               cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
          }
        }

        ++targetCellIndex;
      }
    }
  }

  const bool writeNeedToUpdate = debugMode;
  if (writeNeedToUpdate)
  {
    std::vector<int> needToUpdate(cells.size(), 0);
   
    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
        needToUpdate[cellIndex] = useHalfStepSolution ? 1 : -1;
      }
    }


    MeshVtkWriter< Space, IntTypeCellInfo > meshWriter;

    std::ostringstream needToUpdateFileName;
    needToUpdateFileName << "out/GetData[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(needToUpdateFileName.str(), nodes, cells, needToUpdate);
  }
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::SetCurrCoords(Scalar time, const Scalar* oldCoords)
{
  this->time = time;

  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
  {
    for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
      {
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] = 
          oldCoords[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
      }
    }
  }
  std::fill(inBuffer.begin(), inBuffer.end(), 0);
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState)
{
  const bool writeNeedToUpdate = debugMode;
  if (writeNeedToUpdate)
  {
    std::vector<int> needToUpdate(cells.size(), 0);
    std::vector<int> auxCells(cells.size(), 0);

    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        needToUpdate[cellIndex] = useHalfStepSolution ? 1 : -1;
        auxCells[cellIndex] = auxCell ? 1 : -1;
      }
    }

    MeshVtkWriter< Space, IntTypeCellInfo > meshWriter;

    std::ostringstream needToUpdateFileName;
    needToUpdateFileName << "out/NeedToUpdate[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(needToUpdateFileName.str(), nodes, cells, needToUpdate);

    std::ostringstream auxCellsFileName;
    auxCellsFileName << "out/AuxCell[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(auxCellsFileName.str(), nodes, cells, auxCells);
  }

  this->time = time;

  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;
    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            (useHalfStepSolution ? 
              halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
              cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]) = 
              newCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
          }
        }
        ++targetCellIndex;
      }
    }
  }
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::SetCurrCoords(Scalar time, 
  const Scalar* newCoords, const Scalar* oldCoords, const SolverState& solverState)
{
  const bool writeNeedToUpdate = debugMode;
  if (writeNeedToUpdate)
  {
    std::vector<int> needToUpdate(cells.size(), 0);
    std::vector<int> auxCells(cells.size(), 0);

    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        needToUpdate[cellIndex] = useHalfStepSolution ? 1 : -1;
        auxCells[cellIndex] = auxCell ? 1 : -1;
      }
    }

    MeshVtkWriter< Space, IntTypeCellInfo > meshWriter;

    std::ostringstream needToUpdateFileName;
    needToUpdateFileName << "out/NeedToUpdate[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(needToUpdateFileName.str(), nodes, cells, needToUpdate);

    std::ostringstream auxCellsFileName;
    auxCellsFileName << "out/AuxCell[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(auxCellsFileName.str(), nodes, cells, auxCells);
  }

  this->time = time;

  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;
    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);

        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            if (auxCell || solverState.hierarchyPhase == 1)
            {
              (useHalfStepSolution ? 
                halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
                cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]) =
                  oldCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
            } else 
            {
              // !aux && hierarchyPhase == 0
              halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
                newCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
            }

            if (!auxCell && GetHierarchyLevelsCount() > 1 && solverState.hierarchyPhase == 1)
            {
              bufferCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
                newCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
              inBuffer[cellIndex] = 1;
            }
          }
        }
        ++targetCellIndex;
      }
    }
  }

  // load buffer cell solutions
  if (GetHierarchyLevelsCount() > 1 && solverState.IsPreInitial())
  {
    #pragma omp parallel for
    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      if (inBuffer[cellIndex])
      {
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
              bufferCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];
          }
        }
        inBuffer[cellIndex] = 0;
      }
    }
  }
}

template <typename Space, typename FunctionSpace, typename System>
typename Space::Scalar VolumeMeshCommon<Space, FunctionSpace, System>::
  GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState& solverState, const Scalar* mults)
{
  Scalar maxError = Scalar(0.0);
  const bool writeErrors = debugMode;
  std::vector<Scalar> errors;
  if (writeErrors)
  {
    errors.resize(cells.size(), 0);
  }

  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    for(int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        if (!auxCell && IsCellRegular(cellIndex))
        {
          Scalar cellError = Scalar(0.0);
          for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
          {
            Scalar mult = (mults) ? mults[valueIndex] : Scalar(1.0);
            for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
            {
              cellError += mult * fabs(
                (coords0[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] -
                 coords1[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex])
                * sqrt(cellVolumeIntegrals(functionIndex, functionIndex)));
            }
          }
          if (writeErrors) errors[cellIndex] = cellError;
          #pragma omp flush(maxError)
          if (cellError > maxError)
          {
            #pragma omp critical
            {
              if(cellError > maxError) maxError = cellError;
            }
          }
        }
        ++targetCellIndex;
      }
    }
  }

  /*
    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++)
    {
      Scalar cellError = Scalar(0.0);
      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        Scalar mult = (mults) ? mults[valueIndex] : Scalar(1.0);
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          cellError += mult * fabs(
            (coords0[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] -
            coords1[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex])
            * sqrt(cellVolumeIntegrals[functionsCount * functionIndex + functionIndex]));
        }
      }
      if (cellError > maxError) maxError = cellError;
    } */

  if (writeErrors)
  {
    MeshVtkWriter< Space, ScalarTypeCellInfo > meshWriter;

    std::ostringstream errorsFileName;
    errorsFileName << "out/Errors[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(errorsFileName.str(), nodes, cells, errors);
  }

  return maxError;
}

template <typename Space, typename FunctionSpace, typename System>
typename Space::Scalar VolumeMeshCommon<Space, FunctionSpace, System>::GetTimeStepPrediction()
{
  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
  {
    cellMaxWaveSpeeds[cellIndex] = system.GetMaxWaveSpeed(cellMediumParameters[cellIndex]);
  }

  Scalar minTimeStep = std::numeric_limits<Scalar>::max();
  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
    {
      Scalar cellTimeStep = GetMinHeight(cellIndex) / cellMaxWaveSpeeds[cellIndex];
      minTimeStep = std::min(minTimeStep, cellTimeStep);
    }
  }
  return minTimeStep / Scalar(FunctionSpace::order + 0.5);
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::RebuildTimeHierarchyLevels(IndexType globalStepIndex, 
  bool allowCollisions, bool externalInitialization)
{
  // hierarchy phase count equals 2
  int threadsCount = 1;
  #pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
    {
      threadsCount = omp_get_num_threads();
    }
  }

  std::fill(threadWorkloads.begin(),     threadWorkloads.end(),     0);
  std::fill(threadCellsCount.begin(),    threadCellsCount.end(),    0);
  std::fill(threadCellOffsets.begin(),   threadCellOffsets.end(),   0);
  std::fill(threadSegmentBegins.begin(), threadSegmentBegins.end(), 0);
  std::fill(threadSegmentEnds.begin(),   threadSegmentEnds.end(),   0);

  // build time hierachy levels
  timeHierarchyLevelsManager.Initialize(cells.size(), GetHierarchyLevelsCount(), GetSolverPhasesCount(), externalInitialization);


  Scalar minTimeStep = GetTimeStepPrediction();

  if (globalStepIndex % 100 == 0)
    std::cout << "Min timestep: " << minTimeStep << std::endl;

  if (GetHierarchyLevelsCount() > 1)
  {
    timeHierarchyLevelsManager.BuildTimeHierarchyLevels(this, cellMaxWaveSpeeds, allowCollisions, minTimeStep);

    const bool writeTimeHierarchyLevels = true;
    if (writeTimeHierarchyLevels)
    {
      timeHierarchyLevelsManager.SaveToVtk(this, globalStepIndex);
    }

    std::cout << "Computational cost of " << GetHierarchyLevelsCount() << "-level time hierarchy equals " << 
      timeHierarchyLevelsManager.GetComputationalTotalCost() << std::endl;
  }

  // precomputing
  SolverState solverState(GetHierarchyLevelsCount(), GetSolverPhasesCount());
  for (int globalStepIndex = 0; globalStepIndex < GetMaxHierarchyLevel(); ++globalStepIndex)
  {
    for (int hierarchyPhase = 0; hierarchyPhase < 2; ++hierarchyPhase)
    {
      for (int hierarchyLevel = 0; hierarchyLevel < GetHierarchyLevelsCount(); ++hierarchyLevel)
      {
        solverState.SetState(globalStepIndex, hierarchyPhase, hierarchyLevel);
        IndexType totalCellCount = 0;

        #pragma omp parallel for reduction(+:totalCellCount)
        for (int threadIndex = 0; threadIndex < threadsCount; ++threadIndex)
        {
          IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
          IndexType segmentBegin = threadSegmentBegins[stateIndex] = (threadIndex + 0) * cells.size() / threadsCount;
          IndexType segmentEnd   = threadSegmentEnds[stateIndex]   = (threadIndex + 1) * cells.size() / threadsCount;
          threadCellsCount[stateIndex] = 0;
          for (IndexType cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
          {
            // ***
            if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState))
            {
              ++threadCellsCount[stateIndex];
              if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
              {
                ++threadWorkloads[stateIndex];
              }
            }
          }
          totalCellCount += threadCellsCount[stateIndex];
        }
        hierarchyDimentionsCount[solverState.Index()] = dimsCount * functionsCount * totalCellCount;
        
        // cells per thread balancing
        bool balanced;
        do {
          balanced = false;
          for (int threadIndex = 0; threadIndex + 1 < threadsCount; ++threadIndex)
          {
            IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
            IndexType nextThreadStateIndex = (threadIndex + 1) * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
            while (threadWorkloads[stateIndex] > threadWorkloads[nextThreadStateIndex] + 1)
            {
              --threadSegmentEnds[stateIndex];
              --threadSegmentBegins[nextThreadStateIndex];

              IndexType cellIndex = threadSegmentEnds[stateIndex];
              if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState))
              {
                --threadCellsCount[stateIndex];
                ++threadCellsCount[nextThreadStateIndex];

                if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
                {
                  --threadWorkloads[stateIndex];
                  ++threadWorkloads[nextThreadStateIndex];
                  balanced = true;
                }
              }
            }
            while (threadWorkloads[stateIndex] + 1 < threadWorkloads[nextThreadStateIndex])
            {
              IndexType cellIndex = threadSegmentEnds[stateIndex];
              if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState))
              {
                ++threadCellsCount[stateIndex];
                --threadCellsCount[nextThreadStateIndex];

                if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
                {
                  ++threadWorkloads[stateIndex];
                  --threadWorkloads[nextThreadStateIndex];
                  balanced = true;
                }
              }
              ++threadSegmentEnds[stateIndex];
              ++threadSegmentBegins[nextThreadStateIndex];
            }
          }
        } while (balanced); 
      }
    }
  } 

  for (int globalStepIndex = 0; globalStepIndex < GetMaxHierarchyLevel(); ++globalStepIndex)
  {
    for (int hierarchyPhase = 0; hierarchyPhase < 2; ++hierarchyPhase)
    {
      for (int hierarchyLevel = 0; hierarchyLevel < GetHierarchyLevelsCount(); ++hierarchyLevel)
      {
        solverState.SetState(globalStepIndex, hierarchyPhase, hierarchyLevel);
        for (int threadIndex = 0; threadIndex + 1 < threadsCount; ++threadIndex)
        {
          IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
          IndexType nextThreadStateIndex = (threadIndex + 1) * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
          threadCellOffsets[nextThreadStateIndex] = threadCellOffsets[stateIndex] + threadCellsCount[stateIndex];
        }
      }
    }
  }
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::
  TransformCellSolution(IndexType cellIndex,
    IndexType associatedPermutation[Space::NodesPerCell],
    CellSolution* cellSolution)
{
  bool ordered = true;
  for (IndexType nodeNumber = 0; nodeNumber + 1 < Space::NodesPerCell; ++nodeNumber)
  {
    if (associatedPermutation[nodeNumber] > associatedPermutation[nodeNumber + 1]) ordered = false;
  }
  assert(ordered);
  // TODO for network communication
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::Initialize()
{
  time = Scalar(0.0);
  collisionWidth = Scalar(0.0);
  debugMode = false;

  cellSolutions.resize(cells.size());
  if (GetMaxHierarchyLevel() > 1)
  {
    halfStepCellSolutions.resize(cells.size());
    bufferCellSolutions.resize(cells.size());
  }

  hierarchyDimentionsCount.resize(GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);

  cellMaxWaveSpeeds.resize(cells.size());

  IndexType threadsCount = 1;
  #pragma omp parallel
  {
    threadsCount = omp_get_num_threads();
  }

  isCellAvailable.resize(cells.size(), true);
  inBuffer.resize(cells.size(), 0);
  threadWorkloads.resize(threadsCount     * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadCellsCount.resize(threadsCount    * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadCellOffsets.resize(threadsCount   * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadSegmentBegins.resize(threadsCount * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadSegmentEnds.resize(threadsCount   * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
}

template<typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetCellAverageSolution(IndexType cellIndex) const
{
  typename System::ValueType result;
  std::fill_n(result.values, dimsCount, Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];

      result.values[valueIndex] += basisFunctionCoefficient * cellVolumeAverageIntegrals[functionIndex] / ::Cell<Space>::GetRefCellVolume();
    }
  }
  return result;
}

template<typename Space, typename FunctionSpace, typename System>
bool VolumeMeshCommon<Space, FunctionSpace, System>::IsCellHasDynamicBoundary(IndexType cellIndex) const
{
  bool isCellHasDynamicBoundary = false;

  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
  {
    IndexType correspondingCellIndex = GetCorrespondingCellIndex(cellIndex, faceNumber);
    IndexType correspondingFaceNumber = GetCorrespondingFaceNumber(cellIndex, faceNumber);
    IndexType interactionType = GetInteractionType(cellIndex, faceNumber);

    if (correspondingCellIndex == IndexType(-1) && correspondingFaceNumber == IndexType(-1) && interactionType != IndexType(-1))
    {
      IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);
      if (dynamicContactType != IndexType(-1))
      {
        isCellHasDynamicBoundary = true;
        break;
      }
    }
  }
  return isCellHasDynamicBoundary;
}

template<typename Space, typename FunctionSpace, typename System>
bool VolumeMeshCommon<Space, FunctionSpace, System>::IsReadyForCollisionCell(IndexType cellIndex) const
{
  bool isReadyForCollisionCell = IsCellHasDynamicBoundary(cellIndex);
  
  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
  {
    IndexType correspondingCellIndex = GetCorrespondingCellIndex(cellIndex, faceNumber);

    if (correspondingCellIndex != IndexType(-1) && IsCellHasDynamicBoundary(correspondingCellIndex))
    {
      isReadyForCollisionCell = true;
      break;
    }
  }
  return isReadyForCollisionCell;
}

template<typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::BuildAABBTree(const Vector& boxPoint1, const Vector& boxPoint2)
{
  printf("Building AABB tree for boundary cells\n");
  treeNodeCellIndices.resize(cells.size(), IndexType(-1));

  AABB dynamicContactBox(boxPoint1, boxPoint2);

  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    // if (IsReadyForCollisionCell(cellIndex))
    Vector cellVertices[Space::NodesPerCell];
    GetCellVertices(cellIndex, cellVertices);
    bool needToAdd = false;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (dynamicContactBox.Includes(cellVertices[nodeNumber]))
      {
        needToAdd = true;
        break;
      }
    }

    if (needToAdd)
    {
      treeNodeCellIndices[cellIndex] = aabbTree.InsertNode(GetCellAABB(cellIndex));
      aabbTree.SetUserData(treeNodeCellIndices[cellIndex], cellIndex);
    }
  }

  if (allowDynamicCollisions)
  {
    collisionsInfo.Initialize(cells.size());
  }

  UpdateCollisionsInfo();
}

template<typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::UpdateCollisionsInfo()
{
  if (!allowDynamicCollisions) return;

  collisionsInfo.Clear();

  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    if (!isCellAvailable[cellIndex] || cellMediumParameters[cellIndex].fixed) continue;

    bool isCellContacted = false;
    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
    {
      if (GetCorrespondingCellIndex(cellIndex, faceNumber) == IndexType(-1) &&
        GetCorrespondingFaceNumber(cellIndex, faceNumber) == IndexType(-1))
      {
        IndexType interactionType = GetInteractionType(cellIndex, faceNumber);
        if (interactionType != IndexType(-1))
        {
          IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);
          if (dynamicContactType != IndexType(-1))
          {
            isCellContacted = true;
            break;
          }
        }
      }
    }

    if (isCellContacted)
    {
      ContactFinder< VolumeMeshCommon<Space, FunctionSpace, System> > contactFinder(this, cellIndex);

      const IndexType MaxCollidedCellsCount = 999;
      IndexType collidedCellsCount = 0;
      IndexType collidedCellsPool[MaxCollidedCellsCount];
      contactFinder.Find(collidedCellsPool, &collidedCellsCount, MaxCollidedCellsCount); //second arguments has a pointer

      assert(collidedCellsCount <= MaxCollidedCellsCount);

      if (collidedCellsCount + collisionsInfo.poolSize > collisionsInfo.pool.size())
      {
        collisionsInfo.pool.resize(collidedCellsCount + collisionsInfo.poolSize);
      }

      collisionsInfo.collisionNodes[cellIndex].count = collidedCellsCount;
      collisionsInfo.collisionNodes[cellIndex].offset = collisionsInfo.poolSize;
      collisionsInfo.poolSize += collidedCellsCount;

      if (collidedCellsCount > 0)
      {
        std::copy(collidedCellsPool,
          collidedCellsPool + collidedCellsCount,
          collisionsInfo.pool.data() + collisionsInfo.collisionNodes[cellIndex].offset);
      }
    }
  }
}


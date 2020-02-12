template<typename Space, typename FunctionSpace>
typename Space::IndexType 
  DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::GetSyncDataSize(IndexType dstDomainIndex) const
{
  return syncDataSizes[dstDomainIndex];
}

template<typename Space, typename FunctionSpace>
void DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::
  RebuildTimeHierarchyLevels(IndexType globalStepIndex, bool allowCollisions)
{
  volumeMesh.timeHierarchyLevelsManager.Initialize(
    volumeMesh.cells.size(),
    volumeMesh.GetHierarchyLevelsCount(),
    volumeMesh.GetSolverPhasesCount(), false);

  for (IndexType dstDomainIndex = 0; dstDomainIndex < domainsCount; ++dstDomainIndex)
  {
    for (IndexType cellNumber = 0; cellNumber < transitionInfos[dstDomainIndex].cells.size(); ++cellNumber)
    {
      IndexType cellIndex = volumeMesh.GetCellIndex(transitionInfos[dstDomainIndex].cells[cellNumber].incidentNodes);
      volumeMesh.timeHierarchyLevelsManager.SetLevel(cellIndex, 0);
    }
  }
  volumeMesh.RebuildTimeHierarchyLevels(globalStepIndex, allowCollisions, true);
}

template<typename Space, typename FunctionSpace>
bool DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::
  IsSyncDataEmpty(IndexType dstDomainIndex) const
{
  return transitionInfos[dstDomainIndex].nodesIndices.empty() ||
    transitionInfos[dstDomainIndex].cells.empty() ||
    transitionInfos[dstDomainIndex].transitionNodes.empty();
}

template<typename Space, typename FunctionSpace>
void DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::SetTransitionInfo(
  const std::vector< TransitionInfo<Space> >& transitionInfos)
{
  this->transitionInfos = transitionInfos;
  ComputeSyncDataSizes();

  nodesDictionaries.resize(domainsCount);
  for (IndexType dstDomainIndex = 0; dstDomainIndex < domainsCount; ++dstDomainIndex)
  {
    for (IndexType transitionNodeIndex = 0; transitionNodeIndex < transitionInfos[dstDomainIndex].transitionNodes.size(); ++transitionNodeIndex)
    {
      const typename TransitionInfo<Space>::TransitionNode& transitionNode = transitionInfos[dstDomainIndex].transitionNodes[transitionNodeIndex];
      nodesDictionaries[dstDomainIndex][transitionNode.nativeIndex] = transitionNode.targetIndex;
    }
  }
}

template<typename Space, typename FunctionSpace>
typename DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::BytesHandled
  DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::UpdateDomainData(const char* const data)
{
  IndexType offset = 0;
  // cells
  IndexType cellsCount = *((IndexType*)data);
  offset += sizeof(IndexType);
  UpdateCellsData(cellsCount, (CellSyncData*)(data + offset));
  offset += sizeof(CellSyncData) * cellsCount;

  // nodes
  IndexType nodesCount = *((IndexType*)(data + offset));
  offset += sizeof(IndexType);
  UpdateNodesData(nodesCount, (NodeSyncData*)(data + offset));
  offset += sizeof(NodeSyncData) * nodesCount;
  
  return offset;
}

template<typename Space, typename FunctionSpace>
typename DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::BytesHandled
  DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::BuildSyncData(IndexType dstDomainIndex, char* const data)
{
  *((IndexType*)data) = dstDomainIndex;
  IndexType offset = sizeof(IndexType);

  BuildCellsSyncData(dstDomainIndex, data + offset);
  offset += sizeof(CellSyncData) * transitionInfos[dstDomainIndex].cells.size() + sizeof(IndexType);

  if (sendNodesInfo)
  {
    BuildNodesSyncData(dstDomainIndex, data + offset);
    offset += sizeof(NodeSyncData) * transitionInfos[dstDomainIndex].nodesIndices.size();
  } else
  {
    *((IndexType*)(data + offset)) = 0;
  }
  offset += sizeof(IndexType);

  return offset;
}


template<typename Space, typename FunctionSpace>
void DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::UpdateCellsData(IndexType cellsCount, CellSyncData* const cellsData)
{
  #pragma omp parallel for
  for (int cellNumber = 0; cellNumber < (int)cellsCount; ++cellNumber)
  {
    CellSolution& cellSolution = cellsData[cellNumber].cellSolution;
    IndexType cellIndex = volumeMesh.GetCellIndex(cellsData[cellNumber].cell.incidentNodes);
    assert(cellIndex != IndexType(-1));

    IndexType associatedPermutation[Space::NodesPerCell];
    volumeMesh.GetCellsPairOrientation(cellsData[cellNumber].cell.incidentNodes,
      volumeMesh.cells[cellIndex].incidentNodes, associatedPermutation);

    // ? TODO
    volumeMesh.TransformCellSolution(cellIndex, associatedPermutation, &cellSolution);
    volumeMesh.cellSolutions[cellIndex] = cellSolution;
  }
}

template<typename Space, typename FunctionSpace>
void DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::UpdateNodesData(IndexType nodesCount, const NodeSyncData* const nodesData)
{
  #pragma omp parallel for
  for (int nodeNumber = 0; nodeNumber < (int)nodesCount; ++nodeNumber)
  {
    const Vector& pos = nodesData[nodeNumber].pos;
    IndexType nodeIndex = nodesData[nodeNumber].nodeIndex;
    volumeMesh.nodes[nodeIndex] = Node(pos);
  }
}

template<typename Space, typename FunctionSpace>
void DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::BuildCellsSyncData(IndexType dstDomainIndex, char* const data)
{
  *((IndexType*)data) = (IndexType)transitionInfos[dstDomainIndex].cells.size();
  CellSyncData* cellsData = (CellSyncData*)(data + sizeof(IndexType));

  #pragma omp parallel for
  for (int cellNumber = 0; cellNumber < (int)transitionInfos[dstDomainIndex].cells.size(); ++cellNumber)
  {
    IndexType cellIndex = volumeMesh.GetCellIndex(
      transitionInfos[dstDomainIndex].cells[cellNumber].incidentNodes);
    assert(cellIndex != IndexType(-1));
    CellSyncData& cellSyncData = cellsData[cellNumber];
    cellSyncData.cellSolution = volumeMesh.cellSolutions[cellIndex];
    
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      cellSyncData.cell.incidentNodes[nodeNumber] = nodesDictionaries[dstDomainIndex][volumeMesh.cells[cellIndex].incidentNodes[nodeNumber]];
    }

    // TODO
    // cellSyncData.destroy = 
  }
}

template<typename Space, typename FunctionSpace>
void DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::BuildNodesSyncData(IndexType dstDomainIndex, char* const data)
{
  *((IndexType*)data) = transitionInfos[dstDomainIndex].nodesIndices.size();
  NodeSyncData* nodesData = (NodeSyncData*)(data + sizeof(IndexType));

  #pragma omp parallel for
  for (int nodeNumber = 0; nodeNumber < (int)transitionInfos[dstDomainIndex].nodesIndices.size(); ++nodeNumber)
  {
    IndexType nodeIndex = transitionInfos[dstDomainIndex].nodesIndices[nodeNumber];
    NodeSyncData& nodeSyncData = nodesData[nodeNumber];
    nodeSyncData.pos = volumeMesh.nodes[nodeIndex].pos;
    nodeSyncData.nodeIndex = nodesDictionaries[dstDomainIndex][nodeIndex];
  }
}

template<typename Space, typename FunctionSpace>
void DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::ComputeSyncDataSizes()
{
  for (IndexType dstDomainIndex = 0; dstDomainIndex < domainsCount; ++dstDomainIndex)
  {
    syncDataSizes[dstDomainIndex] += sizeof(IndexType) + // dstDomainIndes
      sizeof(IndexType) + // cells count
      transitionInfos[dstDomainIndex].cells.size() * sizeof(CellSyncData) +
      sizeof(IndexType) + // nodes count
      (sendNodesInfo ? transitionInfos[dstDomainIndex].nodesIndices.size() * sizeof(NodeSyncData) : 0);
  }
}


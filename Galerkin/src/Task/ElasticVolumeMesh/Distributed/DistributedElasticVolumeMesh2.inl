template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::UpdateDomainData(const char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>::UpdateDomainData(data);
  // edges
  IndexType edgesCount = *((IndexType*)(data + offset));
  offset += sizeof(IndexType);
  UpdateEdgesData(edgesCount, (EdgeSyncData*)(data + offset));
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::BuildSyncData(IndexType dstDomainIndex, char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>::BuildSyncData(dstDomainIndex, data);
  if (sendEdgesInfo)
  {
    BuildEdgesSyncData(dstDomainIndex, data + offset);
  } else
  {
    ((IndexType*)(data + offset))[0] = 0;
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::BuildEdgesSyncData(IndexType dstDomainIndex, char* const data)
{
  *((IndexType*)data) = transitionInfos[dstDomainIndex].cells.size() * Space::EdgesPerCell;
  EdgeSyncData* edgesData = (EdgeSyncData*)(data + sizeof(IndexType));

  for (IndexType cellNumber = 0; cellNumber < transitionInfos[dstDomainIndex].cells.size(); ++cellNumber)
  {
    IndexType cellIndex = volumeMesh.GetCellIndex(transitionInfos[dstDomainIndex].cells[cellNumber].incidentNodes);
    assert(cellIndex != IndexType(-1));

    for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; ++edgeNumber)
    {
      EdgeSyncData& edgeSyncData = edgesData[cellNumber * Space::EdgesPerCell + edgeNumber];
      edgeSyncData.interactionType = volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;
      volumeMesh.GetCellEdgeNodes(cellIndex, edgeNumber, edgeSyncData.edge.incidentNodes);
    }
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::UpdateEdgesData(IndexType edgesCount, const EdgeSyncData* const edgesData)
{
  for (IndexType edgeNumber = 0; edgeNumber < edgesCount; ++edgeNumber)
  {
    IndexType interactionType = edgesData[edgeNumber].interactionType;
    // TODO
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::ComputeSyncDataSizes()
{
  DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>::ComputeSyncDataSizes();
  for (IndexType dstDomainIndex = 0; dstDomainIndex < domainsCount; ++dstDomainIndex)
  {
    syncDataSizes[dstDomainIndex] += 
      sizeof(IndexType) + // edges count
      (sendEdgesInfo ? transitionInfos[dstDomainIndex].cells.size() * Space2::EdgesPerCell * sizeof(EdgeSyncData) : 0);
  }
}

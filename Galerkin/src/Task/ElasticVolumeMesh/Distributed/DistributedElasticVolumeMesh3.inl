template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::UpdateDomainData(const char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>::UpdateDomainData(data);
  // faces
  IndexType facesCount = *((IndexType*)(data + offset));
  offset += sizeof(IndexType);
  UpdateFacesData(facesCount, (FaceSyncData*)(data + offset));
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::BuildSyncData(IndexType dstDomainIndex, char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>::BuildSyncData(dstDomainIndex, data);
  if (sendFacesInfo)
  {
    BuildFacesSyncData(dstDomainIndex, data + offset);
  } else
  {
    ((IndexType*)(data + offset))[0] = 0;
  }
}


template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::BuildFacesSyncData(IndexType dstDomainIndex, char* const data)
{
  *((IndexType*)data) = transitionInfos[dstDomainIndex].cells.size() * Space::FacesPerCell;
  FaceSyncData* facesData = (FaceSyncData*)(data + sizeof(IndexType));

  for (IndexType cellNumber = 0; cellNumber < transitionInfos[dstDomainIndex].cells.size(); ++cellNumber)
  {
    IndexType cellIndex = volumeMesh.GetCellIndex(transitionInfos[dstDomainIndex].cells[cellNumber].incidentNodes);
    assert(cellIndex != IndexType(-1));

    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
    {
      FaceSyncData& faceSyncData = facesData[cellNumber * Space::FacesPerCell + faceNumber];
      faceSyncData.interactionType = volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;
      volumeMesh.GetCellFaceNodes(cellIndex, faceNumber, faceSyncData.face.incidentNodes);
    }
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::UpdateFacesData(IndexType facesCount, const FaceSyncData* const facesData)
{
  //assert(0); // TODO
  for (IndexType faceNumber = 0; faceNumber < facesCount; ++faceNumber)
  {
    IndexType interactionType = facesData[faceNumber].interactionType;
    // TODO
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::ComputeSyncDataSizes()
{
  DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>::ComputeSyncDataSizes();
  for (IndexType dstDomainIndex = 0; dstDomainIndex < domainsCount; ++dstDomainIndex)
  {
    syncDataSizes[dstDomainIndex] += 
      sizeof(IndexType) + // faces count
      (sendFacesInfo ? transitionInfos[dstDomainIndex].cells.size() * Space3::FacesPerCell * sizeof(FaceSyncData) : 0);
  }
}


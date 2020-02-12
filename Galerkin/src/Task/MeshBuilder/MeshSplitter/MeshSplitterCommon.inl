template <typename Space>
void MeshSplitterCommon<Space>::Split(const std::vector<IndexType>& cellsDomainIds, 
                                IndexType domainsCount, 
                                std::vector< DistributedMeshIO<Space> >* const domains,
                                bool allowMovement)
{
  this->cellsDomainIds = RenumerateCellsDomainIds(cellsDomainIds);

  domains->resize(domainsCount, DistributedMeshIO<Space>(domainsCount));

  domainsInfos.resize(domainsCount, DomainInfo(domainsCount,
                                               mesh->contactTypesCount, 
                                               mesh->boundaryTypesCount));
  printf("  Building nodes dictionary...\n");
  BuildNodesDictionary(domainsCount, allowMovement);
  SetCorrectSizes(domains);

  printf("  Building nodes...\n");
  // vertices
  for (typename GlobalToLocalIndicesDictionary::const_iterator it = nodesDictionary.begin();
       it != nodesDictionary.end(); ++it)
  {
    (*domains)[it->first.domainIndex].vertices[it->second] = mesh->vertices[it->first.globalIndex];
  }

  printf("  Building cells...\n");
  // cells
  for (IndexType domainIndex = 0; domainIndex < domains->size(); ++domainIndex)
  {
    for (IndexType cellNumber = 0; cellNumber < domainsInfos[domainIndex].cellsIndices.globalIndices.size(); ++cellNumber)
    {
      IndexType cellIndex = domainsInfos[domainIndex].cellsIndices.globalIndices[cellNumber];
      AddCellToDomain(domainIndex, cellIndex, domains);
    }
  }

  printf("  Building detectors...\n");
  std::vector<bool> pointInDomain(domains->size());
  detectorsLocations.resize(mesh->detectorsPositions.size());

  // detectors
  for (IndexType detectorIndex = 0; detectorIndex < mesh->detectorsPositions.size(); ++detectorIndex)
  {
    std::fill(pointInDomain.begin(), pointInDomain.end(), false);
    for (IndexType cellIndex = 0; cellIndex < geomMesh->cells.size(); ++cellIndex)
    {
      IndexType domainIndex = this->cellsDomainIds[cellIndex];
      Vector points[Space::NodesPerCell];
      geomMesh->GetCellVertices(cellIndex, points);
      if (PointInCell<Scalar>(points, mesh->detectorsPositions[detectorIndex]) && !pointInDomain[domainIndex])
      {
        detectorsLocations[detectorIndex] = Location(domainIndex, (*domains)[domainIndex].detectorsPositions.size());
        (*domains)[domainIndex].detectorsPositions.push_back(mesh->detectorsPositions[detectorIndex]);
        pointInDomain[domainIndex] = true;
      }
    }
  }

  printf("  Building corresponding domain pool...\n");
  BuildCellCorrespondingDomainPool(domainsCount);
  printf("  Adding boundaries to domains...\n");
  AddBoundariesToDomains(domainsCount, domains);
  printf("  Adding contacts to domains...\n");
  AddContactsToDomains(domainsCount, domains);

  // sending/receiveing data
  printf("  Building transition cells...\n");
  BuildTransitionCells(domains);
  printf("  Building transition nodes...\n");
  BuildTransitionNodes(domains);
}

template <typename Space>
void MeshSplitterCommon<Space>::BuildNodesDictionary(IndexType domainsCount, bool allowMovement)
{
  for (IndexType cellIndex = 0; cellIndex < geomMesh->cells.size(); ++cellIndex)
  {
    IndexType domainIndex = cellsDomainIds[cellIndex];

    std::vector<IndexType> cellsToHandle;
    cellsToHandle.push_back(cellIndex);

    if (allowMovement)
    {
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
      {
        IndexType nodeIndex = geomMesh->cells[cellIndex].incidentNodes[nodeNumber];
        for (IndexType incidentCellNumber = 0; incidentCellNumber < geomMesh->GetIncidentCellsCount(nodeIndex); ++incidentCellNumber)
        {
          IndexType incidentCellIndex = geomMesh->GetIncidentCellIndex(nodeIndex, incidentCellNumber);

          // auxiliary cell
          if (domainIndex != cellsDomainIds[incidentCellIndex])
          {
            cellsToHandle.push_back(incidentCellIndex);
          }
        }
      }
    }

    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
    {
      IndexType neighbourCellIndex = geomMesh->GetCorrespondingCellIndex(cellIndex, faceNumber);

      if (neighbourCellIndex != IndexType(-1) && domainIndex != cellsDomainIds[neighbourCellIndex])
      {
        cellsToHandle.push_back(neighbourCellIndex);
      }
    }

    for (IndexType cellToHandleNumber = 0; cellToHandleNumber < cellsToHandle.size(); ++cellToHandleNumber)
    {
      domainsInfos[domainIndex].cellsIndices.globalIndices.push_back(cellsToHandle[cellToHandleNumber]);
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
      {
        IndexType nodeIndex = geomMesh->cells[cellsToHandle[cellToHandleNumber]].incidentNodes[nodeNumber];
        domainsInfos[domainIndex].nodesIndices.globalIndices.push_back(nodeIndex);
      }
    }
  }

  for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
  {
    RemoveDuplicates(&domainsInfos[domainIndex].nodesIndices.globalIndices);
    RemoveDuplicates(&domainsInfos[domainIndex].cellsIndices.globalIndices);
  }

  for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
  {
    for (IndexType localNodeIndex = 0; localNodeIndex < domainsInfos[domainIndex].nodesIndices.globalIndices.size(); ++localNodeIndex)
    {
      nodesDictionary[Location(domainIndex, domainsInfos[domainIndex].nodesIndices.globalIndices[localNodeIndex])] = localNodeIndex;
    }
  }
}

template <typename Space>
void MeshSplitterCommon<Space>::SetCorrectSizes(std::vector< DistributedMeshIO<Space> >* const domains)
{
  for (IndexType domainIndex = 0; domainIndex < domainsInfos.size(); ++domainIndex)
  {
    (*domains)[domainIndex].index = domainIndex;
    (*domains)[domainIndex].vertices.resize(domainsInfos[domainIndex].nodesIndices.globalIndices.size());
    (*domains)[domainIndex].indices.reserve(domainsInfos[domainIndex].cellsIndices.globalIndices.size() * Space::NodesPerCell);
  }
}

template <typename Space>
void MeshSplitterCommon<Space>::BuildTransitionCells(std::vector< DistributedMeshIO<Space> >* const domains)
{
  transitionCellDomainInfos.resize(domains->size());
  for (IndexType domainIndex = 0; domainIndex < domains->size(); ++domainIndex)
  {
    transitionCellDomainInfos[domainIndex].transitionCells.resize(domains->size());
  }

  for (IndexType domainIndex = 0; domainIndex < domains->size(); ++domainIndex)
  {
    for (IndexType cellNumber = 0; cellNumber < domainsInfos[domainIndex].cellsIndices.globalIndices.size(); ++cellNumber)
    {
      IndexType cellIndex = domainsInfos[domainIndex].cellsIndices.globalIndices[cellNumber];

      if (domainIndex != cellsDomainIds[cellIndex])
      {
        IndexType neighbourDomainIndex = cellsDomainIds[cellIndex];
        transitionCellDomainInfos[neighbourDomainIndex].transitionCells[domainIndex].cellsIndices.push_back(
          cellsDictionary[Location(neighbourDomainIndex, cellIndex)]);
      }
    }
  }
}

template <typename Space>
void MeshSplitterCommon<Space>::BuildTransitionNodes(std::vector< DistributedMeshIO<Space> >* const domains)
{
  for (IndexType srcDomainIndex = 0; srcDomainIndex < domains->size(); ++srcDomainIndex)
  {
    for (IndexType dstDomainIndex = 0; dstDomainIndex < domains->size(); ++dstDomainIndex)
    {
      // sending nodes
      for (IndexType cellNumber = 0; cellNumber < transitionCellDomainInfos[srcDomainIndex].transitionCells[dstDomainIndex].cellsIndices.size(); ++cellNumber)
      {
        IndexType localCellIndex = transitionCellDomainInfos[srcDomainIndex].transitionCells[dstDomainIndex].cellsIndices[cellNumber];
        IndexType globalCellIndex = domainsInfos[srcDomainIndex].cellsIndices.globalIndices[localCellIndex];
        Cell cell;
        for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
        {
          IndexType globalNodeIndex = geomMesh->cells[globalCellIndex].incidentNodes[nodeNumber];
          IndexType localNodeIndex = nodesDictionary[Location(srcDomainIndex, globalNodeIndex)];
          (*domains)[srcDomainIndex].transitionInfos[dstDomainIndex].nodesIndices.push_back(localNodeIndex);
          cell.incidentNodes[nodeNumber] = localNodeIndex;
        }
        (*domains)[srcDomainIndex].transitionInfos[dstDomainIndex].cells.push_back(cell);
      }
      RemoveDuplicates< IndexType >(&(*domains)[srcDomainIndex].transitionInfos[dstDomainIndex].nodesIndices);
    }
  }

  for (IndexType srcDomainIndex = 0; srcDomainIndex < domains->size(); ++srcDomainIndex)
  {
    for (IndexType dstDomainIndex = 0; dstDomainIndex < domains->size(); ++dstDomainIndex)
    {
      for (IndexType nodeNumber = 0; nodeNumber < (*domains)[srcDomainIndex].transitionInfos[dstDomainIndex].nodesIndices.size(); ++nodeNumber)
      {
        TransitionNode transitionNode;
        transitionNode.nativeIndex = (*domains)[srcDomainIndex].transitionInfos[dstDomainIndex].nodesIndices[nodeNumber];
        IndexType globalNodeIndex = domainsInfos[srcDomainIndex].nodesIndices.globalIndices[transitionNode.nativeIndex];
        transitionNode.targetIndex = nodesDictionary[Location(dstDomainIndex, globalNodeIndex)];
        (*domains)[srcDomainIndex].transitionInfos[dstDomainIndex].transitionNodes.push_back(transitionNode);
      }
    }
  }
}

template <typename Space>
void MeshSplitterCommon<Space>::AddCellToDomain(IndexType domainIndex, IndexType globalCellIndex, 
  std::vector< DistributedMeshIO<Space> >* const domains)
{
  IndexType cellNodeIndices[Space::NodesPerCell];
  geomMesh->GetFixedCellIndices(globalCellIndex, cellNodeIndices);

  cellsDictionary[Location(domainIndex, globalCellIndex)] = (*domains)[domainIndex].GetCellsCount();

  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    IndexType globalNodeIndex = cellNodeIndices[nodeNumber];
    LocalIndex localNodeIndex = nodesDictionary[Location(domainIndex, globalNodeIndex)];
    (*domains)[domainIndex].indices.push_back(localNodeIndex);
  }
  std::sort((*domains)[domainIndex].indices.end() - Space::NodesPerCell, (*domains)[domainIndex].indices.end());
}

template <typename Space>
void MeshSplitterCommon<Space>::BuildCellCorrespondingDomainPool(IndexType domainsCount)
{
  IndexType dataPoolSize = 0;
  cellInfos.resize(mesh->GetCellsCount());

  for (IndexType phase = 0; phase < 2; ++phase)
  {
    if (phase == 1)
    {
      cellsCorrespondingDomainsPool.resize(dataPoolSize);
      IndexType offset = 0;
      for (IndexType cellIndex = 0; cellIndex < mesh->GetCellsCount(); ++cellIndex)
      {
        cellInfos[cellIndex].data = cellInfos[cellIndex].count ? &(cellsCorrespondingDomainsPool[offset]) : 0;
        offset += cellInfos[cellIndex].count;
        cellInfos[cellIndex].count = 0;
      }
    }

    for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
    {
      for (IndexType cellNumber = 0; cellNumber < domainsInfos[domainIndex].cellsIndices.globalIndices.size(); ++cellNumber)
      {
        IndexType cellIndex = domainsInfos[domainIndex].cellsIndices.globalIndices[cellNumber];
        if (phase == 0)
        {
          ++cellInfos[cellIndex].count;
          ++dataPoolSize;
        } else
        {
          cellInfos[cellIndex].data[cellInfos[cellIndex].count++] = domainIndex;
        }
      }
    }
  }
}

template <typename Space>
std::vector<typename Space::IndexType> MeshSplitterCommon<Space>::RenumerateCellsDomainIds(const std::vector<IndexType>& cellsDomainIds)
{
  std::vector<IndexType> updatedDomainIds(cellsDomainIds.size());
  for (IndexType cellIndex = 0; cellIndex < cellsDomainIds.size(); ++cellIndex)
  {
    updatedDomainIds[geomMesh->updatedCellIndices[cellIndex]] = cellsDomainIds[cellIndex];
  }
  return updatedDomainIds;
}

void MeshSplitter<Space2>::AddBoundariesToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains)
{
  typedef std::vector< std::vector<IndexType> > Offsets;
  Offsets boundariesTypeOffsets(domainsCount, std::vector<IndexType>(mesh->boundaryTypesCount));

  for (IndexType phase = 0; phase < 2; ++phase)
  {
    if (phase == 1)
    {
      for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
      {
        if (mesh->boundaryTypesCount > 0)
        {
          std::partial_sum(domainsInfos[domainIndex].boundariesCount.begin(), 
                           domainsInfos[domainIndex].boundariesCount.end(),
                           boundariesTypeOffsets[domainIndex].begin());
        }
        (*domains)[domainIndex].boundaryEdgesCount = domainsInfos[domainIndex].boundariesCount;
        IndexType boundariesCount = std::accumulate(domainsInfos[domainIndex].boundariesCount.begin(),
                                                    domainsInfos[domainIndex].boundariesCount.end()  , 0);
        (*domains)[domainIndex].boundaryEdges.resize(boundariesCount);
        (*domains)[domainIndex].boundaryTypesCount = domainsInfos[domainIndex].boundariesCount.size();
      }
    }
    IndexType offset = 0;

    for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < mesh->boundaryTypesCount; ++boundaryTypeIndex)
    {
      for (IndexType boundaryNumber = 0; boundaryNumber < mesh->boundaryEdgesCount[boundaryTypeIndex]; ++boundaryNumber)
      {
        BoundaryEdge boundaryEdge = mesh->boundaryEdges[offset + boundaryNumber];
        EdgeLocationPair edgeLocationPair = geomMesh->GetEdgeLocation(boundaryEdge);

        for (IndexType pairIndex = 0; pairIndex < 2; ++pairIndex)
        {
          IndexType cellIndex = edgeLocationPair.edges[pairIndex].cellIndex;
          IndexType edgeNumber = edgeLocationPair.edges[pairIndex].edgeNumber;

          if (cellIndex != IndexType(-1) && edgeNumber != IndexType(-1))
          {
            for (IndexType correspondingDomainNumber = 0; correspondingDomainNumber < cellInfos[cellIndex].count; ++correspondingDomainNumber)
            {
              IndexType correspondingDomainIndex = cellInfos[cellIndex].data[correspondingDomainNumber];
              // if it`s not regular cell for this domain, we shouldn`t add any boundaries
              if (correspondingDomainIndex != cellsDomainIds[cellIndex]) continue; 
              if (phase == 0)
              {
                domainsInfos[correspondingDomainIndex].boundariesCount[boundaryTypeIndex]++;
              } else
              {
                BoundaryEdge domainBoundaryEdge;
                for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
                {
                  domainBoundaryEdge.nodeIndices[nodeNumber] = nodesDictionary[Location(correspondingDomainIndex, boundaryEdge.nodeIndices[nodeNumber])];
                }
  
                boundariesTypeOffsets[correspondingDomainIndex][boundaryTypeIndex]--;
                (*domains)[correspondingDomainIndex].boundaryEdges[boundariesTypeOffsets[correspondingDomainIndex][boundaryTypeIndex]] = domainBoundaryEdge;
              }
            }
          }
        }
      }
      offset += mesh->boundaryEdgesCount[boundaryTypeIndex];
    }
  }
}

void MeshSplitter<Space2>::AddContactsToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains)
{
  typedef std::vector< std::vector<IndexType> > Offsets;
  Offsets contactsTypeOffsets(domainsCount, std::vector<IndexType>(mesh->contactTypesCount));

  for (IndexType phase = 0; phase < 2; ++phase)
  {
    if (phase == 1)
    {
      for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
      {
        if (mesh->contactTypesCount > 0)
        {
          std::partial_sum(domainsInfos[domainIndex].contactsCount.begin(), 
                           domainsInfos[domainIndex].contactsCount.end(),
                           contactsTypeOffsets[domainIndex].begin());
        }
        (*domains)[domainIndex].contactEdgesCount = domainsInfos[domainIndex].contactsCount;
        IndexType contactsCount = std::accumulate(domainsInfos[domainIndex].contactsCount.begin(),
                                                  domainsInfos[domainIndex].contactsCount.end()  , 0);
        (*domains)[domainIndex].contactEdges.resize(contactsCount);
        (*domains)[domainIndex].contactTypesCount = domainsInfos[domainIndex].contactsCount.size();
      }
    }
    IndexType offset = 0;

    // contacts exist only if vertices are duplicated
    for (IndexType contactTypeIndex = 0; contactTypeIndex < mesh->contactTypesCount; ++contactTypeIndex)
    {
      for (IndexType contactNumber = 0; contactNumber < mesh->contactEdgesCount[contactTypeIndex]; ++contactNumber)
      {
        EdgePairIndices contactEdge = mesh->contactEdges[offset + contactNumber];
        
        IndexType cellIndices[2];
        for (IndexType pairIndex = 0; pairIndex < 2; ++pairIndex)
        {
          EdgeLocationPair edgeLocationPair = geomMesh->GetEdgeLocation(contactEdge.edges[pairIndex]);
          cellIndices[pairIndex] = IndexType(-1);
          for (IndexType locationIndex = 0; locationIndex < 2; ++locationIndex)
          {
            if (!edgeLocationPair.edges[locationIndex].IsNull())
            {
              cellIndices[pairIndex] = edgeLocationPair.edges[locationIndex].cellIndex;
            }
          }
          assert(cellIndices[pairIndex] != IndexType(-1));
        }

        typedef IndexType CorrespondingDomainNumber;
        for (CorrespondingDomainNumber firstCellDomainNumber = 0; firstCellDomainNumber < cellInfos[cellIndices[0]].count; ++firstCellDomainNumber)
        {
          for (CorrespondingDomainNumber secondCellDomainNumber = 0; secondCellDomainNumber < cellInfos[cellIndices[1]].count; ++secondCellDomainNumber)
          {
            IndexType firstCellDomainIndex = cellInfos[cellIndices[0]].data[firstCellDomainNumber];
            IndexType secondCellDomainIndex = cellInfos[cellIndices[1]].data[secondCellDomainNumber];

            if (firstCellDomainIndex == secondCellDomainIndex)
            {
              if (phase == 0)
              {
                domainsInfos[firstCellDomainIndex].contactsCount[contactTypeIndex]++;
              } else
              if (phase == 1)
              {
                EdgePairIndices edgePair;
                edgePair.edges[0] = EdgeIndices(nodesDictionary[Location(firstCellDomainIndex, contactEdge.edges[0].nodeIndices[0])], 
                                                nodesDictionary[Location(firstCellDomainIndex, contactEdge.edges[0].nodeIndices[1])]);

                edgePair.edges[1] = EdgeIndices(nodesDictionary[Location(firstCellDomainIndex, contactEdge.edges[1].nodeIndices[0])], 
                                                nodesDictionary[Location(firstCellDomainIndex, contactEdge.edges[1].nodeIndices[1])]);

                (*domains)[firstCellDomainIndex].contactEdges[--contactsTypeOffsets[firstCellDomainIndex][contactTypeIndex]] = edgePair;
              }
            }
          }
        }
      }
      offset += mesh->contactEdgesCount[contactTypeIndex];
    }
  }
}

void MeshSplitter<Space3>::AddBoundariesToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains)
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
        (*domains)[domainIndex].boundaryFacesCount = domainsInfos[domainIndex].boundariesCount;
        IndexType boundariesCount = std::accumulate(domainsInfos[domainIndex].boundariesCount.begin(),
          domainsInfos[domainIndex].boundariesCount.end(), 0);
        (*domains)[domainIndex].boundaryFaces.resize(boundariesCount);
        (*domains)[domainIndex].boundaryTypesCount = domainsInfos[domainIndex].boundariesCount.size();
      }
    }
    IndexType offset = 0;

    for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < mesh->boundaryTypesCount; ++boundaryTypeIndex)
    {
      for (IndexType boundaryNumber = 0; boundaryNumber < mesh->boundaryFacesCount[boundaryTypeIndex]; ++boundaryNumber)
      {
        BoundaryFace boundaryFace = mesh->boundaryFaces[offset + boundaryNumber];
        FaceLocationPair faceLocationPair = geomMesh->GetFaceLocation(boundaryFace.nodeIndices);

        for (IndexType pairIndex = 0; pairIndex < 2; ++pairIndex)
        {
          IndexType cellIndex = faceLocationPair.faces[pairIndex].cellIndex;
          IndexType faceNumber = faceLocationPair.faces[pairIndex].faceNumber;

          if (cellIndex != IndexType(-1) && faceNumber != IndexType(-1))
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
                BoundaryFace domainBoundaryFace;
                for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
                {
                  domainBoundaryFace.nodeIndices[nodeNumber] = nodesDictionary[Location(correspondingDomainIndex, boundaryFace.nodeIndices[nodeNumber])];
                }
                boundariesTypeOffsets[correspondingDomainIndex][boundaryTypeIndex]--;
                (*domains)[correspondingDomainIndex].boundaryFaces[boundariesTypeOffsets[correspondingDomainIndex][boundaryTypeIndex]] = domainBoundaryFace;
              }
            }
          }
        }
      }
      offset += mesh->boundaryFacesCount[boundaryTypeIndex];
    }
  }
}

void MeshSplitter<Space3>::AddContactsToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains)
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
        (*domains)[domainIndex].contactFacesCount = domainsInfos[domainIndex].contactsCount;
        IndexType contactsCount = std::accumulate(domainsInfos[domainIndex].contactsCount.begin(),
          domainsInfos[domainIndex].contactsCount.end(), 0);
        (*domains)[domainIndex].contactFaces.resize(contactsCount);
        (*domains)[domainIndex].contactTypesCount = domainsInfos[domainIndex].contactsCount.size();
      }
    }
    IndexType offset = 0;

    // contacts exist only if vertices are duplicated
    for (IndexType contactTypeIndex = 0; contactTypeIndex < mesh->contactTypesCount; ++contactTypeIndex)
    {
      for (IndexType contactNumber = 0; contactNumber < mesh->contactFacesCount[contactTypeIndex]; ++contactNumber)
      {
        FacePairIndices contactFace = mesh->contactFaces[offset + contactNumber];

        IndexType cellIndices[2];
        for (IndexType pairIndex = 0; pairIndex < 2; ++pairIndex)
        {
          FaceLocationPair faceLocationPair = geomMesh->GetFaceLocation(contactFace.faces[pairIndex].nodeIndices);
          cellIndices[pairIndex] = IndexType(-1);
          for (IndexType locationIndex = 0; locationIndex < 2; ++locationIndex)
          {
            if (!faceLocationPair.faces[locationIndex].IsNull())
            {
              cellIndices[pairIndex] = faceLocationPair.faces[locationIndex].cellIndex;
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
                FacePairIndices facePair;
                for (IndexType pair = 0; pair < 2; ++pair)
                {
                  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
                  {
                    facePair.faces[pair].nodeIndices[nodeNumber] = 
                      nodesDictionary[Location(firstCellDomainIndex, contactFace.faces[pair].nodeIndices[nodeNumber])];
                  }
                }
                (*domains)[firstCellDomainIndex].contactFaces[--contactsTypeOffsets[firstCellDomainIndex][contactTypeIndex]] = facePair;
                assert(facePair == contactFace);
              }
            }
          }
        }
      }
      offset += mesh->contactFacesCount[contactTypeIndex];
    }
  }
}

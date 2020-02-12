#include "../DisjointSetBuilder.h"

void GeomMesh<Space3>::BuildAdditionalGroupsInfo(IndexType* submeshNodesCount, IndexType submeshesCount,
  FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
  BoundaryFace*    boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  std::vector<FacePairIndices>  refinedContactFaces;
  std::vector<IndexType>        refinedContactFacesCount;

  //refinedContactFaces.reserve(contactFaces.size());
  refinedContactFacesCount.resize(contactTypesCount);

  IndexType rejectedContactFacesCount = 0;
  IndexType currTypeOffset = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    refinedContactFacesCount[contactTypeIndex] = 0;
    for(IndexType contactFaceIndex = 0; contactFaceIndex < contactFacesCount[contactTypeIndex]; contactFaceIndex++)
    {
      FacePairIndices testPair = contactFaces[currTypeOffset + contactFaceIndex];

      std::sort(testPair.faces[0].nodeIndices, testPair.faces[0].nodeIndices + 3);
      std::sort(testPair.faces[1].nodeIndices, testPair.faces[1].nodeIndices + 3);

      Pair<IndexType> facePair;
      facePair.values[0] = GetFaceIndex(testPair.faces[0].nodeIndices[0], testPair.faces[0].nodeIndices[1], testPair.faces[0].nodeIndices[2]);
      facePair.values[1] = GetFaceIndex(testPair.faces[1].nodeIndices[0], testPair.faces[1].nodeIndices[1], testPair.faces[1].nodeIndices[2]);

      if(facePair.values[0] != IndexType(-1) && facePair.values[1] != IndexType(-1))
      {
        refinedContactFaces.push_back(contactFaces[currTypeOffset + contactFaceIndex]);
        refinedContactFacesCount[contactTypeIndex]++;
      }else
      {
        rejectedContactFacesCount++;
      }
    }
    currTypeOffset += contactFacesCount[contactTypeIndex];
  }
  if(rejectedContactFacesCount > 0)
    printf("rejected contact faces: %d\n", static_cast<int>(rejectedContactFacesCount));

  std::vector<BoundaryFace>  refinedBoundaryFaces;
  std::vector<IndexType>     refinedBoundaryFacesCount;

  // refinedBoundaryFaces.reserve(boundaryFaces.size());
  refinedBoundaryFacesCount.resize(boundaryTypesCount);

  IndexType rejectedBoundaryFacesCount = 0;
  currTypeOffset = 0;
  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    refinedBoundaryFacesCount[boundaryTypeIndex] = 0;
    for(IndexType boundaryFaceIndex = 0; boundaryFaceIndex < boundaryFacesCount[boundaryTypeIndex]; boundaryFaceIndex++)
    {
      BoundaryFace testFace = boundaryFaces[currTypeOffset + boundaryFaceIndex];

      std::sort(testFace.nodeIndices, testFace.nodeIndices + 3);

      IndexType faceIndex = GetFaceIndex(testFace.nodeIndices[0], testFace.nodeIndices[1], testFace.nodeIndices[2]);

      if(faceIndex != IndexType(-1))
      {
        refinedBoundaryFaces.push_back(boundaryFaces[currTypeOffset + boundaryFaceIndex]);
        refinedBoundaryFacesCount[boundaryTypeIndex]++;
      } else
      {
        rejectedBoundaryFacesCount++;
      }
    }
    currTypeOffset += boundaryFacesCount[boundaryTypeIndex];
  }
  if (rejectedBoundaryFacesCount > 0)
    printf("rejected boundary faces: %d\n", static_cast<int>(rejectedBoundaryFacesCount));

  FacePairIndices* contactPtr  = nullptr;
  FaceIndices*     boundaryPtr = nullptr;
  if (refinedContactFaces.size() > 0) contactPtr = refinedContactFaces.data();
  if (refinedBoundaryFaces.size() > 0) boundaryPtr = refinedBoundaryFaces.data();

  BuildNodeContactGroups  (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);
  BuildEdgeContactGroups  (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);
  BuildFaceContactGroups  (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);

  BuildContactNodePairInfos  (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);
  BuildContactEdgePairInfos  (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);
  BuildContactFacePairInfos  (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);

  BuildSmoothContactNormals  (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);
  BuildHardContactNormals    (contactPtr, refinedContactFacesCount.data(), contactTypesCount,  
    boundaryPtr, refinedBoundaryFacesCount.data(), boundaryTypesCount);

  ExpandContactPairInfos();

  BuildSubmeshInfos(submeshNodesCount, submeshesCount);
}


void GeomMesh<Space3>::
  BuildSubmeshInfos(IndexType* submeshNodesCount, IndexType submeshesCount)
{
  if(submeshesCount == 0 || submeshNodesCount == 0)
  {
    this->submeshesCount = 1;
    this->submeshInfos = new SubmeshInfo[1];
    submeshInfos[0].firstNodeIndex = 0;
    submeshInfos[0].nodesCount = GetNodesCount();

    submeshInfos[0].firstEdgeIndex = 0;
    submeshInfos[0].edgesCount = 0;

    submeshInfos[0].firstFaceIndex = 0;
    submeshInfos[0].facesCount = 0;

    submeshInfos[0].firstFaceIndex = 0;
    submeshInfos[0].facesCount = 0;

    submeshInfos[0].firstCellIndex = 0;
    submeshInfos[0].cellsCount = 0;
  } else
  {
    this->submeshesCount = submeshesCount;
    this->submeshInfos = new SubmeshInfo[submeshesCount];

    IndexType currOffset = 0;
    for(IndexType submeshIndex = 0; submeshIndex < submeshesCount; submeshIndex++)
    {
      submeshInfos[submeshIndex].firstNodeIndex = currOffset;
      submeshInfos[submeshIndex].nodesCount = submeshNodesCount[submeshIndex];

      submeshInfos[submeshIndex].firstEdgeIndex = 0;
      submeshInfos[submeshIndex].edgesCount = 0;

      submeshInfos[submeshIndex].firstFaceIndex = 0;
      submeshInfos[submeshIndex].facesCount = 0;

      submeshInfos[submeshIndex].firstFaceIndex = 0;
      submeshInfos[submeshIndex].facesCount = 0;

      submeshInfos[submeshIndex].firstCellIndex = 0;
      submeshInfos[submeshIndex].cellsCount = 0;
      currOffset += submeshNodesCount[submeshIndex];
    }
  }

  IndexType currSubmeshIndex;

  currSubmeshIndex = -1;
  IndexType edgesCount = IndexType(edges.size());
  for(IndexType edgeIndex = 0; (edgeIndex < edgesCount); edgeIndex++)
  {
    IndexType newSubmeshIndex = GetEdgeSubmeshIndex(edgeIndex);

    if (newSubmeshIndex == IndexType(-1))
    {
      std::cout << "Error matching edge submesh indices" << std::endl;
      break;
    }

    if(currSubmeshIndex != newSubmeshIndex)
    {
      submeshInfos[newSubmeshIndex].firstEdgeIndex = edgeIndex;
      submeshInfos[newSubmeshIndex].edgesCount = 0;
      currSubmeshIndex = newSubmeshIndex;
    }

    submeshInfos[currSubmeshIndex].edgesCount++;
  }

  currSubmeshIndex = -1;
  IndexType facesCount = IndexType(faces.size());
  for(IndexType faceIndex = 0; (faceIndex < facesCount); faceIndex++)
  {
    IndexType newSubmeshIndex = GetFaceSubmeshIndex(faceIndex);

    if(newSubmeshIndex == IndexType(-1))
    {
      std::cout << "Error matching face submesh indices" << std::endl;
      break;
    }

    if(currSubmeshIndex != newSubmeshIndex)
    {
      submeshInfos[newSubmeshIndex].firstFaceIndex = faceIndex;
      submeshInfos[newSubmeshIndex].facesCount = 0;
      currSubmeshIndex = newSubmeshIndex;
    }

    submeshInfos[currSubmeshIndex].facesCount++;
  }

  currSubmeshIndex = -1;
  IndexType cellsCount = IndexType(cells.size());
  for(IndexType cellIndex = 0; (cellIndex < cellsCount); cellIndex++)
  {
    IndexType newSubmeshIndex = GetCellSubmeshIndex(cellIndex);

    if(newSubmeshIndex == IndexType(-1))
    {
      std::cout << "Error matching cell submesh indices" << std::endl;
      break;
    }

    if(currSubmeshIndex != newSubmeshIndex)
    {
      submeshInfos[newSubmeshIndex].firstCellIndex = cellIndex;
      submeshInfos[newSubmeshIndex].cellsCount = 0;
      currSubmeshIndex = newSubmeshIndex;
    }

    submeshInfos[currSubmeshIndex].cellsCount++;
  }
}

Space3::IndexType GeomMesh<Space3>::
  GetNodeSubmeshIndex(IndexType nodeIndex)
{
  for(IndexType submeshIndex = 0; submeshIndex < submeshesCount; submeshIndex++)
  {
    if((nodeIndex >= submeshInfos[submeshIndex].firstNodeIndex) &&
       (nodeIndex < submeshInfos[submeshIndex].firstNodeIndex + submeshInfos[submeshIndex].nodesCount))
    {
      return submeshIndex;
    }
  }
  assert(0);
  return -1;
}

Space3::IndexType GeomMesh<Space3>::
  GetEdgeSubmeshIndex(IndexType edgeIndex)
{
  return GetNodeSubmeshIndex(edges[edgeIndex].incidentNodes[0]);
}

Space3::IndexType GeomMesh<Space3>::
  GetFaceSubmeshIndex(IndexType faceIndex)
{
  return GetNodeSubmeshIndex(faces[faceIndex].incidentNodes[0]);
}

Space3::IndexType GeomMesh<Space3>::
  GetCellSubmeshIndex(IndexType cellIndex)
{
  return GetNodeSubmeshIndex(cells[cellIndex].incidentNodes[0]);
}

template<int permutationSize, typename T1, typename T2> 
void TransposeAssociatedPermutation(T1 *referencePermutation, T1 *currentPermutation, T2 *associatedPermutation)
{
  T2 resultPermutation[permutationSize];
  for(int i = 0; i < permutationSize; i++)
  {
    int similarIndex = -1;
    for(int j = 0; j < permutationSize; j++)
    {
      if(currentPermutation[i] == referencePermutation[j]) similarIndex = j;
    }
    assert(similarIndex != (-1));
    resultPermutation[similarIndex] = associatedPermutation[i];
  }
  for(int i = 0; i < permutationSize; i++)
  {
    associatedPermutation[i] = resultPermutation[i];
  }
}

void GeomMesh<Space3>::
  BuildNodeContactGroups(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                         FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  std::vector< Pair<IndexType> > contactNodePairs;
  contactNodePairs.reserve(totalContactFacesCount * 3 + totalBoundaryFacesCount * 3);

  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      Pair<IndexType> newbie;
      newbie.values[0] = contactFaces[i].faces[0].nodeIndices[nodeNumber];
      newbie.values[1] = contactFaces[i].faces[1].nodeIndices[nodeNumber];

      contactNodePairs.push_back(newbie);
    }
  }

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      Pair<IndexType> newbie;
      newbie.values[0] = boundaryFaces[i].nodeIndices[nodeNumber];
      newbie.values[1] = newbie.values[0];

      contactNodePairs.push_back(newbie);
    }
  }
  if(contactNodePairs.size() == 0) return;

  // IndexType nodePairsCount = contactNodePairs.size();

  GroupBuilder<IndexType, IndexType>    contactNodeGroupBuilder;

  contactNodeGroupBuilder.LoadPairs(&(contactNodePairs[0]), contactNodePairs.size());

  contactNodeGroupsCount = contactNodeGroupBuilder.GetGroupsCount();
  contactNodeGroupInfos.resize(contactNodeGroupsCount);

  IndexType nodePairInfosCount = 0;

  IndexType currOffset = 0;
  for(IndexType groupIndex = 0; groupIndex < contactNodeGroupsCount; groupIndex++)
  {
    contactNodeGroupInfos[groupIndex].groupInfoOffset = currOffset;
    contactNodeGroupInfos[groupIndex].pairInfoOffset = nodePairInfosCount;

    IndexType contactPairsCount = contactNodeGroupBuilder.GetGroupSize(groupIndex);
    //contactNodeGroupData.reserve(currOffset + contactPairsCount * 2);

    IndexType nodesFound = 0;
    for(IndexType pairIndex = 0; pairIndex < contactPairsCount; pairIndex++)
    {
      IndexType similarNodeNumber[2];
      IndexType contactIndex = contactNodeGroupBuilder.GetGroupElement(groupIndex, pairIndex);

      for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
      {
        similarNodeNumber[sideIndex] = IndexType(-1);
        IndexType newNodeIndex = contactNodePairs[contactIndex].values[sideIndex];

        for(IndexType nodeNumber = 0; nodeNumber < nodesFound; nodeNumber++)
        {
          IndexType existingNodeIndex = contactNodeGroupData[currOffset + nodeNumber].nodeIndex;
          if(newNodeIndex == existingNodeIndex)
          {
            similarNodeNumber[sideIndex] = nodeNumber;
          }
        }
      }

      if((similarNodeNumber[0] == IndexType(-1) && similarNodeNumber[1] == IndexType(-1))) //both contact sides should be added
      {
        for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
        {
          ContactNode newContactNode;
          newContactNode.nodeIndex = contactNodePairs[contactIndex].values[sideIndex];
          if(newContactNode.nodeIndex == IndexType(-1))
          {
            // int pp = -1;
          }
          //newContactNode.normal = nodes[newContactNode.nodeIndex].pos; //normal is temporary stored in position
          contactNodeGroupData.push_back(newContactNode);
          nodesFound++;

          if(contactNodePairs[contactIndex].values[0] == contactNodePairs[contactIndex].values[1]) break; //boundary node
        }
      }else //add only one side
      {
        for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
        {
          if((similarNodeNumber[sideIndex] == IndexType(-1)) && (similarNodeNumber[(sideIndex + 1) % 2] != IndexType(-1)))
          {
            ContactNode newContactNode;
            newContactNode.nodeIndex = contactNodePairs[contactIndex].values[sideIndex];
            //newContactNode.normal = nodes[newContactNode.nodeIndex].pos; //normal is temporary stored in position

            contactNodeGroupData.push_back(newContactNode);
            nodesFound++;
          }
        }
      }
    }
    currOffset += nodesFound;
    contactNodeGroupInfos[groupIndex].size = nodesFound;
    
    nodePairInfosCount += (contactNodeGroupInfos[groupIndex].size + 1) * (contactNodeGroupInfos[groupIndex].size + 1);
  }

  contactNodePairInfos.resize(nodePairInfosCount);

  for(IndexType nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++)
  {
    nodeContactGroupIndices[nodeIndex] = IndexType(-1);
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex = 0; sideIndex < GetContactNodeGroupSize(groupIndex); sideIndex++)
    {
      nodeContactGroupIndices[GetContactNodeInGroup(groupIndex, sideIndex).nodeIndex] = groupIndex;
    }
  }
}

void GeomMesh<Space3>::
  BuildEdgeContactGroups(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                         FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  std::vector< Pair<IndexType> > contactEdgePairs;
  contactEdgePairs.reserve(totalContactFacesCount * 3 + totalBoundaryFacesCount * 3);

  std::vector<EdgePairIndices > edgePairIndices;
  edgePairIndices.reserve(totalContactFacesCount * 3 + totalBoundaryFacesCount * 3);

  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      EdgePairIndices newbie;
      newbie.edges[0].nodeIndices[0] = contactFaces[i].faces[0].nodeIndices[edgeNumber];
      newbie.edges[0].nodeIndices[1] = contactFaces[i].faces[0].nodeIndices[(edgeNumber + 1) % 3];

      newbie.edges[1].nodeIndices[0] = contactFaces[i].faces[1].nodeIndices[edgeNumber];
      newbie.edges[1].nodeIndices[1] = contactFaces[i].faces[1].nodeIndices[(edgeNumber + 1) % 3];

      edgePairIndices.push_back(newbie);
    }
  }

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      EdgePairIndices newbie;
      newbie.edges[0].nodeIndices[0] = boundaryFaces[i].nodeIndices[edgeNumber];
      newbie.edges[0].nodeIndices[1] = boundaryFaces[i].nodeIndices[(edgeNumber + 1) % 3];
      newbie.edges[1] = newbie.edges[0];

      edgePairIndices.push_back(newbie);
    }
  }

  if (edgePairIndices.size() == 0) return;
  IndexType edgePairsCount = edgePairIndices.size();

  for(IndexType edgePairIndex = 0; edgePairIndex < edgePairsCount; edgePairIndex++)
  {
    EdgePairIndices newbie = edgePairIndices[edgePairIndex];

    std::sort(newbie.edges[0].nodeIndices, newbie.edges[0].nodeIndices + 2);
    std::sort(newbie.edges[1].nodeIndices, newbie.edges[1].nodeIndices + 2);

    Pair<IndexType> edgePair;
    edgePair.values[0] = GetEdgeIndex(newbie.edges[0].nodeIndices[0], newbie.edges[0].nodeIndices[1]);
    edgePair.values[1] = GetEdgeIndex(newbie.edges[1].nodeIndices[0], newbie.edges[1].nodeIndices[1]);

    if(edgePair.values[0] != IndexType(-1) && edgePair.values[1] != IndexType(-1))
    {
      contactEdgePairs.push_back(edgePair);
    } else
    {
      printf("bad contact edge info\n");
    }
  }

  GroupBuilder<IndexType, IndexType>    contactEdgeGroupBuilder;

  contactEdgeGroupBuilder.LoadPairs(&(contactEdgePairs[0]), contactEdgePairs.size());

  contactEdgeGroupsCount = contactEdgeGroupBuilder.GetGroupsCount();
  contactEdgeGroupInfos.resize(contactEdgeGroupsCount);
  IndexType currOffset = 0;
  IndexType edgePairInfosCount = 0;
  for (IndexType groupIndex = 0; groupIndex < contactEdgeGroupsCount; groupIndex++)
  {
    contactEdgeGroupInfos[groupIndex].groupInfoOffset = currOffset;
    contactEdgeGroupInfos[groupIndex].pairInfoOffset = edgePairInfosCount;

    IndexType contactPairsCount = contactEdgeGroupBuilder.GetGroupSize(groupIndex);
    //contactEdgeGroupData.reserve(currOffset + contactPairsCount * 2);

    IndexType edgesFound = 0;
    for(IndexType pairIndex = 0; pairIndex < contactPairsCount; pairIndex++)
    {
      IndexType similarEdgeNumber[2];
      IndexType contactIndex = contactEdgeGroupBuilder.GetGroupElement(groupIndex, pairIndex);

      for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
      {
        similarEdgeNumber[sideIndex] = IndexType(-1);
        IndexType newEdgeIndex = contactEdgePairs[contactIndex].values[sideIndex];

        for(IndexType edgeNumber = 0; edgeNumber < edgesFound; edgeNumber++)
        {
          IndexType existingEdgeIndex = contactEdgeGroupData[currOffset + edgeNumber].edgeIndex;
          if(newEdgeIndex == existingEdgeIndex)
          {
            similarEdgeNumber[sideIndex] = edgeNumber;
          }
        }
      }

      if((similarEdgeNumber[0] == IndexType(-1) && similarEdgeNumber[1] == IndexType(-1))) //both contact sides should be added
      {
        for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
        {
          ContactEdge newContactEdge;
          newContactEdge.edgeIndex = contactEdgePairs[contactIndex].values[sideIndex];
          for(IndexType nodeIndex = 0; nodeIndex < 2; nodeIndex++)
          {
            newContactEdge.incidentNodes[nodeIndex] = edgePairIndices[contactIndex].edges[sideIndex].nodeIndices[nodeIndex];
            //newContactEdge.normals[nodeIndex] = nodes[newContactEdge.incidentNodes[nodeIndex]].pos; //normal is temporary stored in position
          }
          contactEdgeGroupData.push_back(newContactEdge);
          edgesFound++;

          if(contactEdgePairs[contactIndex].values[0] == contactEdgePairs[contactIndex].values[1]) break; //boundary edge
        }
      } else //contact needs rotation
      {
        for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
        {
          if((similarEdgeNumber[sideIndex] == IndexType(-1)) && (similarEdgeNumber[(sideIndex + 1) % 2] != IndexType(-1)))
          {
            ContactEdge newContactEdge;
            newContactEdge.edgeIndex = contactEdgePairs[contactIndex].values[sideIndex];

            for(IndexType nodeIndex = 0; nodeIndex < 2; nodeIndex++)
            {
              newContactEdge.incidentNodes[nodeIndex] = edgePairIndices[contactIndex].edges[sideIndex].nodeIndices[nodeIndex];
              //newContactEdge.normals[nodeIndex] = nodes[newContactEdge.incidentNodes[nodeIndex]].pos; //normal is temporary stored in position
            }

            TransposeAssociatedPermutation<2>(
              contactEdgeGroupData[currOffset + similarEdgeNumber[(sideIndex + 1) % 2]].incidentNodes,
              edgePairIndices[contactIndex].edges[(sideIndex + 1) % 2].nodeIndices,
              newContactEdge.incidentNodes);

            /*TransposeAssociatedPermutation<2>(
              contactEdgeGroupData[currOffset + similarEdgeNumber[(sideIndex + 1) % 2]].incidentNodes,
              edgePairIndices[contactIndex].edges[(sideIndex + 1) % 2].nodeIndices,
              newContactEdge.normals);*/

            contactEdgeGroupData.push_back(newContactEdge);
            edgesFound++;
          }
        }
      }
    }
    currOffset += edgesFound;
    contactEdgeGroupInfos[groupIndex].size = edgesFound;

    edgePairInfosCount += (contactEdgeGroupInfos[groupIndex].size + 1) * (contactEdgeGroupInfos[groupIndex].size + 1);
  }

  contactEdgePairInfos.resize(edgePairInfosCount);

  for(IndexType edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++)
  {
    edgeContactGroupIndices[edgeIndex] = IndexType(-1);
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex = 0; sideIndex < GetContactEdgeGroupSize(groupIndex); sideIndex++)
    {
      edgeContactGroupIndices[GetContactEdgeInGroup(groupIndex, sideIndex).edgeIndex] = groupIndex;
    }
  }
}

void GeomMesh<Space3>::
  BuildFaceContactGroups(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                         FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  std::vector< Pair<IndexType> > contactFacePairs;
  contactFacePairs.reserve(totalContactFacesCount);

  std::vector<FacePairIndices> facePairIndices;
  facePairIndices.reserve(totalContactFacesCount + totalBoundaryFacesCount);

  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    FacePairIndices newbie;

    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      newbie.faces[0].nodeIndices[nodeNumber] = contactFaces[i].faces[0].nodeIndices[nodeNumber];
      newbie.faces[1].nodeIndices[nodeNumber] = contactFaces[i].faces[1].nodeIndices[nodeNumber];
    }

    facePairIndices.push_back(newbie);
  }

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    FacePairIndices newbie;

    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      newbie.faces[0].nodeIndices[nodeNumber] = boundaryFaces[i].nodeIndices[nodeNumber];
      newbie.faces[1].nodeIndices[nodeNumber] = boundaryFaces[i].nodeIndices[nodeNumber];
    }

    facePairIndices.push_back(newbie);
  }

  if (facePairIndices.size() == 0) return;
  IndexType facePairsCount = facePairIndices.size();

  for(IndexType facePairIndex = 0; facePairIndex < facePairsCount; facePairIndex++)
  {
    FacePairIndices newbie = facePairIndices[facePairIndex];

    std::sort(newbie.faces[0].nodeIndices, newbie.faces[0].nodeIndices + 3);
    std::sort(newbie.faces[1].nodeIndices, newbie.faces[1].nodeIndices + 3);

    Pair<IndexType> facePair;
    facePair.values[0] = GetFaceIndex(newbie.faces[0].nodeIndices);
    facePair.values[1] = GetFaceIndex(newbie.faces[1].nodeIndices);

    if(facePair.values[0] != IndexType(-1) && facePair.values[1] != IndexType(-1))
    {
      contactFacePairs.push_back(facePair);
    } else
    {
      printf("bad contact face info\n");
    }
  }

  GroupBuilder<IndexType, IndexType>    contactFaceGroupBuilder;

  contactFaceGroupBuilder.LoadPairs(&(contactFacePairs[0]), contactFacePairs.size());

  contactFaceGroupsCount = contactFaceGroupBuilder.GetGroupsCount();
  contactFaceGroupInfos.resize(contactFaceGroupsCount);
  IndexType currOffset = 0;
  IndexType facePairInfosCount = 0;
  for(IndexType groupIndex = 0; groupIndex < contactFaceGroupsCount; groupIndex++)
  {
    contactFaceGroupInfos[groupIndex].groupInfoOffset = currOffset;
    contactFaceGroupInfos[groupIndex].pairInfoOffset = facePairInfosCount;

    IndexType contactPairsCount = contactFaceGroupBuilder.GetGroupSize(groupIndex);
    //contactFaceGroupData.reserve(currOffset + contactPairsCount * 2);

    IndexType facesFound = 0;
    for(IndexType pairIndex = 0; pairIndex < contactPairsCount; pairIndex++)
    {
      IndexType similarFaceNumber[2];
      IndexType contactIndex = contactFaceGroupBuilder.GetGroupElement(groupIndex, pairIndex);

      for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
      {
        similarFaceNumber[sideIndex] = IndexType(-1);
        IndexType newFaceIndex = contactFacePairs[contactIndex].values[sideIndex];

        for(IndexType faceNumber = 0; faceNumber < facesFound; faceNumber++)
        {
          IndexType existingFaceIndex = contactFaceGroupData[currOffset + faceNumber].faceIndex;
          if(newFaceIndex == existingFaceIndex)
          {
            similarFaceNumber[sideIndex] = faceNumber;
          }
        }
      }

      if((similarFaceNumber[0] == IndexType(-1) && similarFaceNumber[1] == IndexType(-1))) //both contact sides should be added
      {
        for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
        {
          ContactFace newContactFace;
          newContactFace.faceIndex = contactFacePairs[contactIndex].values[sideIndex];
          for(IndexType nodeIndex = 0; nodeIndex < 3; nodeIndex++)
          {
            newContactFace.incidentNodes[nodeIndex] = facePairIndices[contactIndex].faces[sideIndex].nodeIndices[nodeIndex];
            //newContactFace.incidentNodeGroups[nodeIndex] = nodes[newContactFace.incidentNodes[nodeIndex]].contactGroupIndex;
            //newContactFace.normals[nodeIndex] = nodes[newContactFace.incidentNodes[nodeIndex]].pos; //normal is temporary stored in position
          }
          contactFaceGroupData.push_back(newContactFace);
          facesFound++;

          if(contactFacePairs[contactIndex].values[0] == contactFacePairs[contactIndex].values[1]) break; //boundary face
        }
      }else
      {
        for(IndexType sideIndex = 0; sideIndex < 2; sideIndex++)
        {
          if((similarFaceNumber[sideIndex] == IndexType(-1)) && (similarFaceNumber[(sideIndex + 1) % 2] != IndexType(-1)))
          {
            ContactFace newContactFace;
            newContactFace.faceIndex = contactFacePairs[contactIndex].values[sideIndex];

            for(IndexType nodeIndex = 0; nodeIndex < 3; nodeIndex++)
            {
              newContactFace.incidentNodes[nodeIndex] = facePairIndices[contactIndex].faces[sideIndex].nodeIndices[nodeIndex];
              //newContactFace.incidentNodeGroups[nodeIndex] = nodes[newContactFace.incidentNodes[nodeIndex]].contactGroupIndex;
              //newContactFace.normals[nodeIndex] = nodes[newContactFace.incidentNodes[nodeIndex]].pos; //normal is temporary stored in position
            }

            TransposeAssociatedPermutation<3>(
              contactFaceGroupData[currOffset + similarFaceNumber[(sideIndex + 1) % 2]].incidentNodes,
              facePairIndices[contactIndex].faces[(sideIndex + 1) % 2].nodeIndices,
              newContactFace.incidentNodes);

            /*TransposeAssociatedPermutation<3>(
              contactFaceGroupData[currOffset + similarFaceNumber[(sideIndex + 1) % 2]].incidentNodes,
              facePairIndices[contactIndex].faces[(sideIndex + 1) % 2].nodeIndices,
              newContactFace.normals);*/

            contactFaceGroupData.push_back(newContactFace);
            facesFound++;
          }
        }
      }
    }
    currOffset += facesFound;
    contactFaceGroupInfos[groupIndex].size = facesFound;

    facePairInfosCount += (contactFaceGroupInfos[groupIndex].size + 1) * (contactFaceGroupInfos[groupIndex].size + 1);
  }

  contactFacePairInfos.resize(facePairInfosCount);

  for(IndexType faceIndex = 0; faceIndex < faces.size(); faceIndex++)
  {
    faceContactGroupIndices[faceIndex] = IndexType(-1);
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex = 0; sideIndex < GetContactFaceGroupSize(groupIndex); sideIndex++)
    {
      faceContactGroupIndices[GetContactFaceInGroup(groupIndex, sideIndex).faceIndex] = groupIndex;
    }
  }
}

template<typename Scalar, typename Vector3>
Vector3 GetWeightedNormal(Vector3 v0, Vector3 v1)
{
  Vector3 normal = v0.GetNorm() ^ v1.GetNorm();
  Scalar sinAng = normal.Len();
  if(sinAng > Scalar(1.0)) sinAng = Scalar(1.0); //precision bugs workaround
  if(sinAng < Scalar(-1.0)) sinAng = Scalar(-1.0);
  Scalar ang = asin(sinAng);
  return normal.GetNorm() * ang;
}

void GeomMesh<Space3>::
  BuildContactNodePairInfos(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                            FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex0 = 0; sideIndex0 < GetContactNodeGroupSize(groupIndex) + 1; sideIndex0++)
    {
      for(IndexType sideIndex1 = 0; sideIndex1 < GetContactNodeGroupSize(groupIndex) + 1; sideIndex1++)
      {
        GetContactNodePairInfo(groupIndex, sideIndex0, sideIndex1).typeIndex = IndexType(-1);
        GetContactNodePairInfo(groupIndex, sideIndex0, sideIndex1).normal = Vector::zeroVector();
      }
    }
  }

  IndexType currContactType = 0;
  IndexType currContactTypeEnd = contactFacesCount[currContactType];
  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    for(;(i >= currContactTypeEnd) && (currContactType < contactTypesCount);
        currContactTypeEnd += contactFacesCount[++currContactType]);

    IndexType groupSideIndices[2][3];
    IndexType groupIndices[2][3];
    for(IndexType faceNumber = 0; faceNumber < 2; faceNumber++)
    {
      for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
      {
        IndexType nodeIndex = contactFaces[i].faces[faceNumber].nodeIndices[nodeNumber];

        IndexType groupIndex = nodeContactGroupIndices[nodeIndex];
        assert(groupIndex != -1);
        groupIndices[faceNumber][nodeNumber] = groupIndex;
        for(IndexType sideIndex = 0; sideIndex < GetContactNodeGroupSize(groupIndex); sideIndex++)
        {
          if(GetContactNodeInGroup(groupIndex, sideIndex).nodeIndex == nodeIndex)
          {
            groupSideIndices[faceNumber][nodeNumber] = sideIndex;
          }
        }
      }
    }

    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      assert(groupIndices[0][nodeNumber] == groupIndices[1][nodeNumber]); //both sides should be in the same group
      IndexType groupIndex = groupIndices[0][nodeNumber];

      if(groupSideIndices[0][nodeNumber] != groupSideIndices[1][nodeNumber])
      {
        ContactNodePairInfo& contactInfo01 = GetContactNodePairInfo(groupIndex, 
          groupSideIndices[0][nodeNumber], groupSideIndices[1][nodeNumber]);
        if((contactInfo01.typeIndex == IndexType(-1)) || 
          (currContactType < contactInfo01.typeIndex && contactInfo01.typeIndex != IndexType(-1)))
        {
          contactInfo01.typeIndex = currContactType;
        }
        ContactNodePairInfo& contactInfo10 = GetContactNodePairInfo(groupIndex, groupSideIndices[1][nodeNumber], 
          groupSideIndices[0][nodeNumber]);
        if((contactInfo10.typeIndex == IndexType(-1)) || 
          (currContactType < contactInfo10.typeIndex && contactInfo10.typeIndex != IndexType(-1)))
        {
          contactInfo10.typeIndex = currContactType;
        }
      }
    }
  }

  IndexType currBoundaryType = 0;
  IndexType currBoundaryTypeEnd = boundaryFacesCount[currBoundaryType];

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    for(;(i >= currBoundaryTypeEnd) && (currBoundaryType < boundaryTypesCount);
        currBoundaryTypeEnd += boundaryFacesCount[++currBoundaryType]);

    IndexType groupSideIndices[3];
    IndexType groupIndices[3];

    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      IndexType nodeIndex = boundaryFaces[i].nodeIndices[nodeNumber];

      groupIndices[nodeNumber] = nodeContactGroupIndices[nodeIndex];
      assert(groupIndices[nodeNumber] != -1);
      for(IndexType sideIndex = 0; sideIndex < GetContactNodeGroupSize(groupIndices[nodeNumber]); sideIndex++)
      {
        if(GetContactNodeInGroup(groupIndices[nodeNumber], sideIndex).nodeIndex == nodeIndex)
        {
          groupSideIndices[nodeNumber] = sideIndex;
        }
      }
    }

    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      ContactNodePairInfo& contactInfo0b = GetContactNodePairInfo(groupIndices[nodeNumber], groupSideIndices[nodeNumber], -1);
      if ((contactInfo0b.typeIndex == IndexType(-1)) || 
        (currBoundaryType < contactInfo0b.typeIndex && contactInfo0b.typeIndex != IndexType(-1)))
      {
        contactInfo0b.typeIndex = currBoundaryType;
      }

      ContactNodePairInfo& contactInfob0 = GetContactNodePairInfo(groupIndices[nodeNumber], -1, groupSideIndices[nodeNumber]);
      if ((contactInfob0.typeIndex == IndexType(-1)) || 
        (currBoundaryType < contactInfob0.typeIndex && contactInfob0.typeIndex != IndexType(-1)))
      {
        contactInfob0.typeIndex = currBoundaryType;
      }
    }
  }
}

void GeomMesh<Space3>::
  BuildContactEdgePairInfos(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                            FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex0 = 0; sideIndex0 < GetContactEdgeGroupSize(groupIndex) + 1; sideIndex0++)
    {
      for(IndexType sideIndex1 = 0; sideIndex1 < GetContactEdgeGroupSize(groupIndex) + 1; sideIndex1++)
      {
        GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1).typeIndex = IndexType(-1);
        GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1).normals[0] = Vector::zeroVector();
        GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1).normals[1] = Vector::zeroVector();
      }
    }
  }

  IndexType currContactType = 0;
  IndexType currContactTypeEnd = contactFacesCount[currContactType];
  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    for(;(i >= currContactTypeEnd) && (currContactType < contactTypesCount);
        currContactTypeEnd += contactFacesCount[++currContactType]);

    IndexType groupSideIndices[2][3];
    IndexType groupIndices[2][3];
    for(IndexType faceNumber = 0; faceNumber < 2; faceNumber++)
    {
      for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
      {
        IndexType nodeIndices[2];
        nodeIndices[0] = contactFaces[i].faces[faceNumber].nodeIndices[(edgeNumber + 0) % 3];
        nodeIndices[1] = contactFaces[i].faces[faceNumber].nodeIndices[(edgeNumber + 1) % 3];

        std::sort(nodeIndices, nodeIndices + 2);

        IndexType edgeIndex = GetEdgeIndex(nodeIndices[0], nodeIndices[1]);
        assert(edgeIndex != IndexType(-1));

        IndexType groupIndex = edgeContactGroupIndices[edgeIndex];
        assert(groupIndex != IndexType(-1));
        groupIndices[faceNumber][edgeNumber] = groupIndex;

        groupSideIndices[faceNumber][edgeNumber] = GetEdgeSideInGroup(groupIndex, edgeIndex);
      }
    }

    for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      assert(groupIndices[0][(edgeNumber + 0) % 3] == groupIndices[1][(edgeNumber + 0) % 3]); //both sides should be in the same group
      assert(groupIndices[0][(edgeNumber + 1) % 3] == groupIndices[1][(edgeNumber + 1) % 3]); //both sides should be in the same group
      assert(groupIndices[0][(edgeNumber + 2) % 3] == groupIndices[1][(edgeNumber + 2) % 3]); //both sides should be in the same group

      IndexType groupIndex = groupIndices[0][edgeNumber];

      if (groupSideIndices[0][edgeNumber] != groupSideIndices[1][edgeNumber])
      {
        ContactEdgePairInfo& contactInfo01 = GetContactEdgePairInfo(groupIndex, 
          groupSideIndices[0][edgeNumber], groupSideIndices[1][edgeNumber]);

        if (contactInfo01.typeIndex == IndexType(-1) || 
          (currContactType < contactInfo01.typeIndex && contactInfo01.typeIndex != IndexType(-1)))
        {
          contactInfo01.typeIndex = currContactType;
        }

        ContactEdgePairInfo& contactInfo10 = GetContactEdgePairInfo(groupIndex, 
          groupSideIndices[1][edgeNumber], groupSideIndices[0][edgeNumber]);

        if (contactInfo10.typeIndex == IndexType(-1) || 
          (currContactType < contactInfo10.typeIndex && contactInfo10.typeIndex != IndexType(-1)))
        {
          contactInfo10.typeIndex = currContactType;
        }
      }
    }
  }

  IndexType currBoundaryType = 0;
  IndexType currBoundaryTypeEnd = boundaryFacesCount[currBoundaryType];

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    for(;(i >= currBoundaryTypeEnd) && (currBoundaryType < boundaryTypesCount);
        currBoundaryTypeEnd += boundaryFacesCount[++currBoundaryType]);

    IndexType groupSideIndices[3];
    IndexType groupIndices[3];

    for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      IndexType nodeIndices[2];
      nodeIndices[0] = boundaryFaces[i].nodeIndices[(edgeNumber + 0) % 3];
      nodeIndices[1] = boundaryFaces[i].nodeIndices[(edgeNumber + 1) % 3];

      std::sort(nodeIndices, nodeIndices + 2);

      IndexType edgeIndex = GetEdgeIndex(nodeIndices[0], nodeIndices[1]);
      assert(edgeIndex != IndexType(-1));

      groupIndices[edgeNumber] = edgeContactGroupIndices[edgeIndex];
      assert(groupIndices[edgeNumber] != IndexType(-1));

      groupSideIndices[edgeNumber] = GetEdgeSideInGroup(groupIndices[edgeNumber], edgeIndex);
    }

    for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      ContactEdgePairInfo& contactInfo0b = GetContactEdgePairInfo(groupIndices[edgeNumber], groupSideIndices[edgeNumber], -1);
      if (contactInfo0b.typeIndex == IndexType(-1) || 
        (currBoundaryType < contactInfo0b.typeIndex && contactInfo0b.typeIndex != IndexType(-1)))
      {
        contactInfo0b.typeIndex = currBoundaryType;
      }

      ContactEdgePairInfo& contactInfob0 = GetContactEdgePairInfo(groupIndices[edgeNumber], -1, groupSideIndices[edgeNumber]);
      if ( contactInfob0.typeIndex == IndexType(-1) || 
        (currBoundaryType < contactInfob0.typeIndex && contactInfob0.typeIndex != IndexType(-1)) )
      {
        contactInfob0.typeIndex = currBoundaryType;
      }
    }
  }
}

void GeomMesh<Space3>::
  BuildContactFacePairInfos(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                            FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount  = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex0 = 0; sideIndex0 < GetContactFaceGroupSize(groupIndex) + 1; sideIndex0++)
    {
      for(IndexType sideIndex1 = 0; sideIndex1 < GetContactFaceGroupSize(groupIndex) + 1; sideIndex1++)
      {
        GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1).typeIndex = IndexType(-1);
        GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1).normals[0] = Vector::zeroVector();
        GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1).normals[1] = Vector::zeroVector();
        GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1).normals[2] = Vector::zeroVector();
      }
    }
  }

  IndexType currContactType = 0;
  IndexType currContactTypeEnd = contactFacesCount[currContactType];
  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    for(;(i >= currContactTypeEnd) && (currContactType < contactTypesCount);
        currContactTypeEnd += contactFacesCount[++currContactType]);

    IndexType groupSideIndex[2];
    IndexType groupIndex[2];
    for(IndexType faceNumber = 0; faceNumber < 2; faceNumber++)
    {
      IndexType nodeIndices[3];
      nodeIndices[0] = contactFaces[i].faces[faceNumber].nodeIndices[0];
      nodeIndices[1] = contactFaces[i].faces[faceNumber].nodeIndices[1];
      nodeIndices[2] = contactFaces[i].faces[faceNumber].nodeIndices[2];

      std::sort(nodeIndices, nodeIndices + 3);

      IndexType faceIndex = GetFaceIndex(nodeIndices[0], nodeIndices[1], nodeIndices[2]);
      assert(faceIndex != IndexType(-1));

      groupIndex[faceNumber] = faceContactGroupIndices[faceIndex];
      assert(groupIndex[faceNumber] != IndexType(-1));

      groupSideIndex[faceNumber] = GetFaceSideInGroup(groupIndex[faceNumber], faceIndex);
    }

    assert(groupIndex[0] == groupIndex[1]); //both sides should be in the same group

    if(groupSideIndex[0] != groupSideIndex[1])
    {
      ContactFacePairInfo& contactInfo01 = GetContactFacePairInfo(groupIndex[0], groupSideIndex[0], groupSideIndex[1]);
      if ( contactInfo01.typeIndex == IndexType(-1) || 
        (currContactType < contactInfo01.typeIndex && contactInfo01.typeIndex != IndexType(-1)) )
      {
        contactInfo01.typeIndex = currContactType;
      }

      ContactFacePairInfo& contactInfo10 = GetContactFacePairInfo(groupIndex[0], groupSideIndex[1], groupSideIndex[0]);
      if ( contactInfo10.typeIndex == IndexType(-1) || 
        (currContactType < contactInfo10.typeIndex && contactInfo10.typeIndex != IndexType(-1)) )
      {
        contactInfo10.typeIndex = currContactType;
      }
    }
  }

  IndexType currBoundaryType = 0;
  IndexType currBoundaryTypeEnd = boundaryFacesCount[currBoundaryType];

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    for(;(i >= currBoundaryTypeEnd) && (currBoundaryType < boundaryTypesCount);
        currBoundaryTypeEnd += boundaryFacesCount[++currBoundaryType]);

    IndexType groupSideIndex;
    IndexType groupIndex;

    IndexType nodeIndices[3];
    nodeIndices[0] = boundaryFaces[i].nodeIndices[0];
    nodeIndices[1] = boundaryFaces[i].nodeIndices[1];
    nodeIndices[2] = boundaryFaces[i].nodeIndices[2];

    std::sort(nodeIndices, nodeIndices + 3);

    IndexType faceIndex = GetFaceIndex(nodeIndices[0], nodeIndices[1], nodeIndices[2]);
    assert(faceIndex != IndexType(-1));

    groupIndex = faceContactGroupIndices[faceIndex];
    assert(groupIndex != IndexType(-1));

    groupSideIndex = GetFaceSideInGroup(groupIndex, faceIndex);
    assert(groupSideIndex != IndexType(-1));

    ContactFacePairInfo& contactInfo0b = GetContactFacePairInfo(groupIndex, groupSideIndex, -1);
    if (contactInfo0b.typeIndex == IndexType(-1) || 
      (currBoundaryType < contactInfo0b.typeIndex && contactInfo0b.typeIndex != IndexType(-1)))
    {
      contactInfo0b.typeIndex = currBoundaryType;
    }

    ContactFacePairInfo& contactInfob0 = GetContactFacePairInfo(groupIndex, -1, groupSideIndex);
    if ( contactInfob0.typeIndex == IndexType(-1) || 
      (currBoundaryType < contactInfob0.typeIndex && contactInfob0.typeIndex != IndexType(-1)) )
    {
      contactInfob0.typeIndex = currBoundaryType;
    }
  }
}

void GeomMesh<Space3>::
  BuildSmoothContactNormals(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                            FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  std::vector<Vector> nodeNormals;
  nodeNormals.resize(GetNodesCount());

  for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
  {
    nodeNormals[nodeIndex] = Vector::zeroVector();
  }

  IndexType currContactType = 0;
  IndexType currContactTypeEnd = contactFacesCount[currContactType];
  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    for(;(i >= currContactTypeEnd) && (currContactType < contactTypesCount);
        currContactTypeEnd += contactFacesCount[++currContactType]);

    FacePairIndices newbie = contactFaces[i];

    for(IndexType faceNumber = 0; faceNumber < 2; faceNumber++)
    {
      std::sort(newbie.faces[faceNumber].nodeIndices, newbie.faces[faceNumber].nodeIndices + 3);
      IndexType faceIndex = GetFaceIndex(
        newbie.faces[faceNumber].nodeIndices[0], 
        newbie.faces[faceNumber].nodeIndices[1], 
        newbie.faces[faceNumber].nodeIndices[2]);
      assert(faceIndex != -1);
      Vector outerNormal = GetFaceOuterNormal(faceIndex, nullptr);

      for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
      {
        IndexType nodeIndex = newbie.faces[faceNumber].nodeIndices[nodeNumber];

        Vector p0 = nodes[newbie.faces[faceNumber].nodeIndices[(nodeNumber + 0) % 3]].pos;
        Vector p1 = nodes[newbie.faces[faceNumber].nodeIndices[(nodeNumber + 1) % 3]].pos;
        Vector p2 = nodes[newbie.faces[faceNumber].nodeIndices[(nodeNumber + 2) % 3]].pos;
        Vector norm = GetWeightedNormal<Scalar, Vector>((p1 - p0), (p2 - p0));
        if (norm * outerNormal < Scalar(0)) norm = -norm;

        nodeNormals[nodeIndex] += norm;
      }
    }

    if(i == currContactTypeEnd - 1)
    {
      /*for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
      {
        Scalar len = absolute(nodeNormals[nodeIndex]);

        if(len != 0)
          nodeNormals[nodeIndex] *= Scalar(1.0) / len;
        else
          nodeNormals[nodeIndex] = Vector3(-10, 0, 0);
      }*/

      //node normals
      for(IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex0 = 0; sideIndex0 < GetContactNodeGroupSize(groupIndex); sideIndex0++)
        {
          for(IndexType sideIndex1 = 0; sideIndex1 < GetContactNodeGroupSize(groupIndex); sideIndex1++)
          {
            ContactNodePairInfo& nodePairInfo = GetContactNodePairInfo(groupIndex, sideIndex0, sideIndex1);

            if(nodePairInfo.typeIndex == currContactType)
            {
              nodePairInfo.normal = nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex0).nodeIndex] - 
                                    nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex1).nodeIndex];
            }
          }
        }
      }

      //edge normals
      for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex0 = 0; sideIndex0 < GetContactEdgeGroupSize(groupIndex); sideIndex0++)
        {
          for(IndexType sideIndex1 = 0; sideIndex1 < GetContactEdgeGroupSize(groupIndex); sideIndex1++)
          {
            ContactEdgePairInfo& edgePairInfo = GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1);

            if(edgePairInfo.typeIndex == currContactType)
            {
              edgePairInfo.normals[0] = 
                nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex0).incidentNodes[0]] - 
                nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex1).incidentNodes[0]];
              edgePairInfo.normals[1] = 
                nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex0).incidentNodes[1]] - 
                nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex1).incidentNodes[1]];
            }
          }
        }
      }

      //face normals
      for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex0 = 0; sideIndex0 < GetContactFaceGroupSize(groupIndex); sideIndex0++) //+1's are border points
        {
          for(IndexType sideIndex1 = 0; sideIndex1 < GetContactFaceGroupSize(groupIndex); sideIndex1++)
          {
            ContactFacePairInfo& facePairInfo = GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1);

            if(facePairInfo.typeIndex == currContactType)
            {
              facePairInfo.normals[0] = 
                nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex0).incidentNodes[0]] - 
                nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex1).incidentNodes[0]];
              facePairInfo.normals[1] = 
                nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex0).incidentNodes[1]] - 
                nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex1).incidentNodes[1]];
              facePairInfo.normals[2] = 
                nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex0).incidentNodes[2]] - 
                nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex1).incidentNodes[2]];
            }
          }
        }
      }

      for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
      {
        nodeNormals[nodeIndex] = Vector::zeroVector();
      }
    }
  }

  for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
  {
    nodeNormals[nodeIndex] = Vector::zeroVector();
  }

  IndexType currBoundaryType = 0;
  IndexType currBoundaryTypeEnd = boundaryFacesCount[currBoundaryType];

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    for(;(i >= currBoundaryTypeEnd) && (currBoundaryType < boundaryTypesCount);
        currBoundaryTypeEnd += boundaryFacesCount[++currBoundaryType]);


    FaceIndices newbie = boundaryFaces[i];

    std::sort(newbie.nodeIndices, newbie.nodeIndices + 3);
    IndexType faceIndex = GetFaceIndex(
      newbie.nodeIndices[0], 
      newbie.nodeIndices[1], 
      newbie.nodeIndices[2]);
    assert(faceIndex != -1);

    Vector outerNormal = GetFaceOuterNormal(faceIndex,  nullptr);

    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      IndexType nodeIndex = boundaryFaces[i].nodeIndices[nodeNumber];
      Vector p0 = nodes[boundaryFaces[i].nodeIndices[(nodeNumber + 0) % 3]].pos;
      Vector p1 = nodes[boundaryFaces[i].nodeIndices[(nodeNumber + 1) % 3]].pos;
      Vector p2 = nodes[boundaryFaces[i].nodeIndices[(nodeNumber + 2) % 3]].pos;
      Vector norm = GetWeightedNormal<Scalar, Vector>((p1 - p0), (p2 - p0));
      if (norm * outerNormal < Scalar(0)) norm = -norm;

      nodeNormals[nodeIndex] += norm;
    }

    if(i == currBoundaryTypeEnd - 1)
    {
      /*for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
      {
        Scalar len = absolute(nodeNormals[nodeIndex]);

        if(len != 0)
          nodeNormals[nodeIndex] *= Scalar(1.0) / len;
        else
          nodeNormals[nodeIndex] = Vector3(-10, 0, 0);
      }*/
      //node boundary normals
      for(IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex = 0; sideIndex < GetContactNodeGroupSize(groupIndex); sideIndex++)
        {
          ContactNodePairInfo& nodePairInfo0b = GetContactNodePairInfo(groupIndex, sideIndex, -1);

          if(nodePairInfo0b.typeIndex == currBoundaryType)
          {
            nodePairInfo0b.normal = nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex).nodeIndex];
          }

          ContactNodePairInfo& nodePairInfob0 = GetContactNodePairInfo(groupIndex, -1, sideIndex);
          if(nodePairInfob0.typeIndex == currBoundaryType)
          {
            nodePairInfob0.normal = -nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex).nodeIndex];
          }
        }
      }

      //edge boundary normals
      for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex = 0; sideIndex < GetContactEdgeGroupSize(groupIndex); sideIndex++)
        {
          ContactEdgePairInfo& edgePairInfo0b = GetContactEdgePairInfo(groupIndex, sideIndex, -1);

          if (edgePairInfo0b.typeIndex == currBoundaryType)
          {
            edgePairInfo0b.normals[0] =  nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex).incidentNodes[0]];
            edgePairInfo0b.normals[1] =  nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex).incidentNodes[1]];
          }

          ContactEdgePairInfo& edgePairInfob0 = GetContactEdgePairInfo(groupIndex, -1, sideIndex);
          if (edgePairInfob0.typeIndex == currBoundaryType)
          {
            edgePairInfob0.normals[0] = -nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex).incidentNodes[0]];
            edgePairInfob0.normals[1] = -nodeNormals[GetContactEdgeInGroup(groupIndex, sideIndex).incidentNodes[1]];
          }
        }
      }

      //face boundary normals
      for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex = 0; sideIndex < GetContactFaceGroupSize(groupIndex); sideIndex++)
        {
          ContactFacePairInfo& facePairInfo0b = GetContactFacePairInfo(groupIndex, sideIndex, -1);

          if(facePairInfo0b.typeIndex == currBoundaryType)
          {
            facePairInfo0b.normals[0] =  nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex).incidentNodes[0]];
            facePairInfo0b.normals[1] =  nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex).incidentNodes[1]];
            facePairInfo0b.normals[2] =  nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex).incidentNodes[2]];
          }

          ContactFacePairInfo& facePairInfob0 = GetContactFacePairInfo(groupIndex, -1, sideIndex);
          if(facePairInfob0.typeIndex == currBoundaryType)
          {
            facePairInfob0.normals[0] = -nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex).incidentNodes[0]];
            facePairInfob0.normals[1] = -nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex).incidentNodes[1]];
            facePairInfob0.normals[2] = -nodeNormals[GetContactFaceInGroup(groupIndex, sideIndex).incidentNodes[2]];
          }
        }
      }

      for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
      {
        nodeNormals[nodeIndex] = Vector::zeroVector();
      }
    }
  }
}

void GeomMesh<Space3>::
  BuildHardContactNormals(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                          FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType totalContactFacesCount = 0;
  IndexType totalBoundaryFacesCount = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    totalContactFacesCount += contactFacesCount[contactTypeIndex];
  }

  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    totalBoundaryFacesCount += boundaryFacesCount[boundaryTypeIndex];
  }

  std::vector<Vector> nodeNormals;
  nodeNormals.resize(GetNodesCount());

  std::vector<Vector> edgeNormals;
  edgeNormals.resize(GetEdgesCount());

  std::vector<Vector> faceNormals;
  faceNormals.resize(GetFacesCount());

  for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
  {
    nodeNormals[nodeIndex] = Vector::zeroVector();
  }

  for(IndexType edgeIndex = 0; edgeIndex < GetEdgesCount(); edgeIndex++)
  {
    edgeNormals[edgeIndex] = Vector::zeroVector();
  }

  for(IndexType faceIndex = 0; faceIndex < GetFacesCount(); faceIndex++)
  {
    faceNormals[faceIndex] = Vector::zeroVector();
  }

  IndexType currContactType = 0;
  IndexType currContactTypeEnd = contactFacesCount[currContactType];
  for(IndexType i = 0; i < totalContactFacesCount; i++)
  {
    for(;(i >= currContactTypeEnd) && (currContactType < contactTypesCount);
        currContactTypeEnd += contactFacesCount[++currContactType]);

    FacePairIndices newbie = contactFaces[i];

    for(IndexType faceNumber = 0; faceNumber < 2; faceNumber++)
    {
      std::sort(newbie.faces[faceNumber].nodeIndices, newbie.faces[faceNumber].nodeIndices + 3);
      IndexType faceIndex = GetFaceIndex(
        newbie.faces[faceNumber].nodeIndices[0], 
        newbie.faces[faceNumber].nodeIndices[1], 
        newbie.faces[faceNumber].nodeIndices[2]);
      assert(faceIndex != -1);

      Vector outerNormal = GetFaceOuterNormal(faceIndex, nullptr);

      for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
      {
        IndexType nodeIndex = newbie.faces[faceNumber].nodeIndices[nodeNumber];

        Vector p0 = nodes[newbie.faces[faceNumber].nodeIndices[(nodeNumber + 0) % 3]].pos;
        Vector p1 = nodes[newbie.faces[faceNumber].nodeIndices[(nodeNumber + 1) % 3]].pos;
        Vector p2 = nodes[newbie.faces[faceNumber].nodeIndices[(nodeNumber + 2) % 3]].pos;
        Vector norm = GetWeightedNormal<Scalar, Vector>((p1 - p0), (p2 - p0));
        if (norm * outerNormal < Scalar(0)) norm = -norm;
        
        nodeNormals[nodeIndex] += norm;
      }

      for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
      {
        IndexType incidentNodes[2];
        incidentNodes[0] = newbie.faces[faceNumber].nodeIndices[(edgeNumber + 0) % 3];
        incidentNodes[1] = newbie.faces[faceNumber].nodeIndices[(edgeNumber + 1) % 3];
        std::sort(incidentNodes, incidentNodes + 2);

        IndexType edgeIndex = GetEdgeIndex(incidentNodes[0], incidentNodes[1]);
        assert(edgeIndex != IndexType(-1));

        edgeNormals[edgeIndex] += outerNormal;
      }

      faceNormals[faceIndex] += outerNormal;
    }

    if(i == currContactTypeEnd - 1)
    {
      //node normals
      for(IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex0 = 0; sideIndex0 < GetContactNodeGroupSize(groupIndex); sideIndex0++)
        {
          for(IndexType sideIndex1 = 0; sideIndex1 < GetContactNodeGroupSize(groupIndex); sideIndex1++)
          {
            ContactNodePairInfo& nodePairInfo = GetContactNodePairInfo(groupIndex, sideIndex0, sideIndex1);

            if(nodePairInfo.typeIndex == currContactType)
            {
              nodePairInfo.normal = nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex0).nodeIndex] - 
                                    nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex1).nodeIndex];
            }
          }
        }
      }

      //edge normals
      for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex0 = 0; sideIndex0 < GetContactEdgeGroupSize(groupIndex); sideIndex0++)
        {
          for(IndexType sideIndex1 = 0; sideIndex1 < GetContactEdgeGroupSize(groupIndex); sideIndex1++)
          {
            ContactEdgePairInfo& edgePairInfo = GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1);

            if(edgePairInfo.typeIndex == currContactType)
            {
              edgePairInfo.normals[0] = 
                edgeNormals[GetContactEdgeInGroup(groupIndex, sideIndex0).edgeIndex] - 
                edgeNormals[GetContactEdgeInGroup(groupIndex, sideIndex1).edgeIndex];
              edgePairInfo.normals[1] = edgePairInfo.normals[0];
            }
          }
        }
      }

      //face normals
      for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex0 = 0; sideIndex0 < GetContactFaceGroupSize(groupIndex); sideIndex0++) //+1's are border points
        {
          for(IndexType sideIndex1 = 0; sideIndex1 < GetContactFaceGroupSize(groupIndex); sideIndex1++)
          {
            ContactFacePairInfo& facePairInfo = GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1);

            if(facePairInfo.typeIndex == currContactType)
            {
              facePairInfo.normals[0] = 
                faceNormals[GetContactFaceInGroup(groupIndex, sideIndex0).faceIndex] - 
                faceNormals[GetContactFaceInGroup(groupIndex, sideIndex1).faceIndex];
              facePairInfo.normals[1] = facePairInfo.normals[0];
              facePairInfo.normals[2] = facePairInfo.normals[0];
            }
          }
        }
      }

      for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
      {
        nodeNormals[nodeIndex] = Vector::zeroVector();
      }
      for(IndexType edgeIndex = 0; edgeIndex < GetEdgesCount(); edgeIndex++)
      {
        edgeNormals[edgeIndex] = Vector::zeroVector();
      }
      for(IndexType faceIndex = 0; faceIndex < GetFacesCount(); faceIndex++)
      {
        faceNormals[faceIndex] = Vector::zeroVector();
      }
    }
  }

  for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
  {
    nodeNormals[nodeIndex] = Vector::zeroVector();
  }
  for(IndexType edgeIndex = 0; edgeIndex < GetEdgesCount(); edgeIndex++)
  {
    edgeNormals[edgeIndex] = Vector::zeroVector();
  }
  for(IndexType faceIndex = 0; faceIndex < GetFacesCount(); faceIndex++)
  {
    faceNormals[faceIndex] = Vector::zeroVector();
  }

  IndexType currBoundaryType = 0;
  IndexType currBoundaryTypeEnd = boundaryFacesCount[currBoundaryType];

  for(IndexType i = 0; i < totalBoundaryFacesCount; i++)
  {
    for(;(i >= currBoundaryTypeEnd) && (currBoundaryType < boundaryTypesCount);
        currBoundaryTypeEnd += boundaryFacesCount[++currBoundaryType]);


    FaceIndices newbie = boundaryFaces[i];

    std::sort(newbie.nodeIndices, newbie.nodeIndices + 3);
    IndexType faceIndex = GetFaceIndex(
      newbie.nodeIndices[0], 
      newbie.nodeIndices[1], 
      newbie.nodeIndices[2]);
    assert(faceIndex != -1);

    Vector outerNormal = GetFaceOuterNormal(faceIndex,  nullptr);

    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      IndexType nodeIndex = boundaryFaces[i].nodeIndices[nodeNumber];
      Vector p0 = nodes[boundaryFaces[i].nodeIndices[(nodeNumber + 0) % 3]].pos;
      Vector p1 = nodes[boundaryFaces[i].nodeIndices[(nodeNumber + 1) % 3]].pos;
      Vector p2 = nodes[boundaryFaces[i].nodeIndices[(nodeNumber + 2) % 3]].pos;
      Vector norm = GetWeightedNormal<Scalar, Vector>((p1 - p0), (p2 - p0));
      if (norm * outerNormal < Scalar(0)) norm = -norm;

      nodeNormals[nodeIndex] += norm;
    }

    for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      IndexType incidentNodes[2];
      incidentNodes[0] = boundaryFaces[i].nodeIndices[(edgeNumber + 0) % 3];
      incidentNodes[1] = boundaryFaces[i].nodeIndices[(edgeNumber + 1) % 3];
      std::sort(incidentNodes, incidentNodes + 2);

      IndexType edgeIndex = GetEdgeIndex(incidentNodes[0], incidentNodes[1]);
      assert(edgeIndex != IndexType(-1));

      edgeNormals[edgeIndex] += outerNormal;
    }

    faceNormals[faceIndex] += outerNormal;

    if(i == currBoundaryTypeEnd - 1)
    {
      //node boundary normals
      for(IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex = 0; sideIndex < GetContactNodeGroupSize(groupIndex); sideIndex++)
        {
          ContactNodePairInfo& nodePairInfo0b = GetContactNodePairInfo(groupIndex, sideIndex, -1);

          if(nodePairInfo0b.typeIndex == currBoundaryType)
          {
            nodePairInfo0b.normal = nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex).nodeIndex];
          }

          ContactNodePairInfo& nodePairInfob0 = GetContactNodePairInfo(groupIndex, -1, sideIndex);
          if(nodePairInfob0.typeIndex == currBoundaryType)
          {
            nodePairInfob0.normal = -nodeNormals[GetContactNodeInGroup(groupIndex, sideIndex).nodeIndex];
          }
        }
      }

      //edge boundary normals
      for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex = 0; sideIndex < GetContactEdgeGroupSize(groupIndex); sideIndex++)
        {
          ContactEdgePairInfo& edgePairInfo0b = GetContactEdgePairInfo(groupIndex, sideIndex, -1);

          if(edgePairInfo0b.typeIndex == currBoundaryType)
          {
            edgePairInfo0b.normals[0] =  edgeNormals[GetContactEdgeInGroup(groupIndex, sideIndex).edgeIndex];
            edgePairInfo0b.normals[1] =  edgePairInfo0b.normals[0];
          }

          ContactEdgePairInfo& edgePairInfob0 = GetContactEdgePairInfo(groupIndex, -1, sideIndex);
          if(edgePairInfob0.typeIndex == currBoundaryType)
          {
            edgePairInfob0.normals[0] = -edgeNormals[GetContactEdgeInGroup(groupIndex, sideIndex).edgeIndex];
            edgePairInfob0.normals[1] = edgePairInfob0.normals[0];
          }
        }
      }

      //face boundary normals
      for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
      {
        for(IndexType sideIndex = 0; sideIndex < GetContactFaceGroupSize(groupIndex); sideIndex++)
        {
          ContactFacePairInfo& facePairInfo0b = GetContactFacePairInfo(groupIndex, sideIndex, -1);

          if(facePairInfo0b.typeIndex == currBoundaryType)
          {
            facePairInfo0b.normals[0] = faceNormals[GetContactFaceInGroup(groupIndex, sideIndex).faceIndex];
            facePairInfo0b.normals[1] = facePairInfo0b.normals[0];
            facePairInfo0b.normals[2] = facePairInfo0b.normals[0];
          }

          ContactFacePairInfo& facePairInfob0 = GetContactFacePairInfo(groupIndex, -1, sideIndex);
          if(facePairInfob0.typeIndex == currBoundaryType)
          {
            facePairInfob0.normals[0] = -faceNormals[GetContactFaceInGroup(groupIndex, sideIndex).faceIndex];
            facePairInfob0.normals[1] = facePairInfob0.normals[0];
            facePairInfob0.normals[2] = facePairInfob0.normals[0];
          }
        }
      }

      for(IndexType nodeIndex = 0; nodeIndex < GetNodesCount(); nodeIndex++)
      {
        nodeNormals[nodeIndex] = Vector::zeroVector();
      }
      for(IndexType edgeIndex = 0; edgeIndex < GetEdgesCount(); edgeIndex++)
      {
        edgeNormals[edgeIndex] = Vector::zeroVector();
      }
      for(IndexType faceIndex = 0; faceIndex < GetFacesCount(); ++faceIndex)
      {
        faceNormals[faceIndex] = Vector::zeroVector();
      }
    }
  }
}

void GeomMesh<Space3>::
  ExpandContactPairInfos()
{
  //expanding contact node matrix data
  for(IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
  {
    for(IndexType midIndex = 0; midIndex < GetContactNodeGroupSize(groupIndex); midIndex++)
    {
      for(IndexType sideIndex0 = 0; sideIndex0 < GetContactNodeGroupSize(groupIndex); sideIndex0++)
      {
        for(IndexType sideIndex1 = 0; sideIndex1 < GetContactNodeGroupSize(groupIndex); sideIndex1++)
        {
          ContactNodePairInfo& straightPath  = GetContactNodePairInfo(groupIndex, sideIndex0, sideIndex1);
          ContactNodePairInfo& splitPath0    = GetContactNodePairInfo(groupIndex, sideIndex0, midIndex);
          ContactNodePairInfo& splitPath1    = GetContactNodePairInfo(groupIndex, midIndex, sideIndex1);

          if (splitPath0.typeIndex == IndexType(-1)) continue;
          if (splitPath1.typeIndex == IndexType(-1)) continue;

          IndexType resultType;
          Vector resultNormal;

          if (splitPath0.typeIndex > splitPath1.typeIndex)
          {
            resultType = splitPath0.typeIndex;
            resultNormal = splitPath0.normal;
          } else
          if (splitPath0.typeIndex < splitPath1.typeIndex)
          {
            resultType = splitPath1.typeIndex;
            resultNormal = splitPath1.normal;
          } else
          {
            assert(splitPath0.typeIndex == splitPath1.typeIndex);

            resultType = splitPath0.typeIndex;
            resultNormal = splitPath0.normal + splitPath1.normal;
          }

          if( (straightPath.typeIndex == IndexType(-1) || straightPath.typeIndex > resultType) && resultType != IndexType(-1) )
          {
            straightPath.typeIndex = resultType;
            straightPath.normal = resultNormal;
          }
        }
      }
    }
  }

  for (IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
  {
    for (IndexType sideIndex = 0; sideIndex < GetContactNodeGroupSize(groupIndex); sideIndex++)
    {
      ContactNodePairInfo& generalInfo = GetContactNodePairInfo(groupIndex, -1, -1);
      ContactNodePairInfo& currentInfo = GetContactNodePairInfo(groupIndex, sideIndex, -1);

      if ( (currentInfo.typeIndex < generalInfo.typeIndex || generalInfo.typeIndex == IndexType(-1)) &&
         currentInfo.typeIndex != IndexType(-1) )
      {
        generalInfo.typeIndex = currentInfo.typeIndex;
        generalInfo.normal = Vector::zeroVector();
      }
      if (generalInfo.typeIndex == currentInfo.typeIndex)
      {
        generalInfo.normal += currentInfo.normal;
      }
    }
  }

  for (IndexType groupIndex = 0; groupIndex < GetContactNodeGroupsCount(); groupIndex++)
  {
    for (IndexType sideIndex0 = 0; sideIndex0 < GetContactNodeGroupSize(groupIndex) + 1; sideIndex0++)
    {
      for (IndexType sideIndex1 = 0; sideIndex1 < GetContactNodeGroupSize(groupIndex) + 1; sideIndex1++)
      {
        ContactNodePairInfo& currentInfo = GetContactNodePairInfo(groupIndex, sideIndex0, sideIndex1);
        if (currentInfo.typeIndex != IndexType(-1))
        {
          Scalar len = (currentInfo.normal).Len();

          if(len != 0)
            currentInfo.normal *= Scalar(1.0) / len;
          else
            currentInfo.normal = Vector::zeroVector();
        }
      }
    }
  }

  //expanding contact edge matrix data
  for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
  {
    for(IndexType midIndex = 0; midIndex < GetContactEdgeGroupSize(groupIndex); midIndex++)
    {
      for(IndexType sideIndex0 = 0; sideIndex0 < GetContactEdgeGroupSize(groupIndex); sideIndex0++)
      {
        for(IndexType sideIndex1 = 0; sideIndex1 < GetContactEdgeGroupSize(groupIndex); sideIndex1++)
        {
          ContactEdgePairInfo& straightPath = GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1);
          ContactEdgePairInfo& splitPath0   = GetContactEdgePairInfo(groupIndex, sideIndex0, midIndex);
          ContactEdgePairInfo& splitPath1   = GetContactEdgePairInfo(groupIndex, midIndex, sideIndex1);

          if(splitPath0.typeIndex == IndexType(-1)) continue;
          if(splitPath1.typeIndex == IndexType(-1)) continue;

          IndexType resultType;
          Vector resultNormals[2];

          if (splitPath0.typeIndex > splitPath1.typeIndex)
          {
            resultType = splitPath0.typeIndex;
            for(IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
              resultNormals[normalNumber] = splitPath0.normals[normalNumber];
          } else
          if (splitPath0.typeIndex < splitPath1.typeIndex)
          {
            resultType = splitPath1.typeIndex;
            for(IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
              resultNormals[normalNumber] = splitPath1.normals[normalNumber];
          } else
          {
            assert(splitPath0.typeIndex == splitPath1.typeIndex);

            resultType = splitPath0.typeIndex;
            for(IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
              resultNormals[normalNumber] = splitPath0.normals[normalNumber] + splitPath1.normals[normalNumber];
          }

          if ( (straightPath.typeIndex == IndexType(-1) || straightPath.typeIndex > resultType) && resultType != IndexType(-1) )
          {
            straightPath.typeIndex = resultType;
            for(IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
              straightPath.normals[normalNumber] = resultNormals[normalNumber];
          }
        }
      }
    }
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex = 0; sideIndex < GetContactEdgeGroupSize(groupIndex); sideIndex++)
    {
      ContactEdgePairInfo& generalInfo = GetContactEdgePairInfo(groupIndex, -1, -1);
      ContactEdgePairInfo& currentInfo = GetContactEdgePairInfo(groupIndex, sideIndex, -1);

      if ( (currentInfo.typeIndex < generalInfo.typeIndex || generalInfo.typeIndex == IndexType(-1)) &&
         currentInfo.typeIndex != IndexType(-1) )
      {
        generalInfo.typeIndex = currentInfo.typeIndex;
        for (IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
          generalInfo.normals[normalNumber] = Vector::zeroVector();
      }
      if (generalInfo.typeIndex == currentInfo.typeIndex)
      {
        for(IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
          generalInfo.normals[normalNumber] += currentInfo.normals[normalNumber];
      }
    }
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactEdgeGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex0 = 0; sideIndex0 < GetContactEdgeGroupSize(groupIndex) + 1; sideIndex0++)
    {
      for(IndexType sideIndex1 = 0; sideIndex1 < GetContactEdgeGroupSize(groupIndex) + 1; sideIndex1++)
      {
        ContactEdgePairInfo& currentInfo = GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1);
        if (currentInfo.typeIndex != IndexType(-1))
        {
          for(IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
          {
            Scalar len = (currentInfo.normals[normalNumber]).Len();

            if(len != 0)
              currentInfo.normals[normalNumber] *= Scalar(1.0) / len;
            else
              currentInfo.normals[normalNumber] = Vector::zeroVector();
          }
        }
      }
    }
  }

  //expanding contact face matrix data
  for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
  {
    for(IndexType midIndex = 0; midIndex < GetContactFaceGroupSize(groupIndex); midIndex++)
    {
      for(IndexType sideIndex0 = 0; sideIndex0 < GetContactFaceGroupSize(groupIndex); sideIndex0++)
      {
        for(IndexType sideIndex1 = 0; sideIndex1 < GetContactFaceGroupSize(groupIndex); sideIndex1++)
        {
          ContactFacePairInfo& straightPath = GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1);
          ContactFacePairInfo& splitPath0   = GetContactFacePairInfo(groupIndex, sideIndex0, midIndex);
          ContactFacePairInfo& splitPath1   = GetContactFacePairInfo(groupIndex, midIndex, sideIndex1);

          if(splitPath0.typeIndex == IndexType(-1)) continue;
          if(splitPath1.typeIndex == IndexType(-1)) continue;

          IndexType resultType;
          Vector resultNormals[3];

          if (splitPath0.typeIndex > splitPath1.typeIndex)
          {
            resultType = splitPath0.typeIndex;
            for(IndexType normalNumber = 0; normalNumber < 3; normalNumber++)
              resultNormals[normalNumber] = splitPath0.normals[normalNumber];
          } else
          if (splitPath0.typeIndex < splitPath1.typeIndex)
          {
            resultType = splitPath1.typeIndex;
            for(IndexType normalNumber = 0; normalNumber < 3; normalNumber++)
              resultNormals[normalNumber] = splitPath1.normals[normalNumber];
          } else
          {
            assert(splitPath0.typeIndex == splitPath1.typeIndex);
            resultType = splitPath0.typeIndex;
            for(IndexType normalNumber = 0; normalNumber < 3; normalNumber++)
              resultNormals[normalNumber] = splitPath0.normals[normalNumber] + splitPath1.normals[normalNumber];
          }
          
          if ( (straightPath.typeIndex == IndexType(-1) || straightPath.typeIndex > resultType) && resultType != IndexType(-1) )
          {
            straightPath.typeIndex = resultType;
            for(IndexType normalNumber = 0; normalNumber < 3; normalNumber++)
              straightPath.normals[normalNumber] = resultNormals[normalNumber];
          }
        }
      }
    }
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex = 0; sideIndex < GetContactFaceGroupSize(groupIndex); sideIndex++)
    {
      ContactFacePairInfo& generalInfo = GetContactFacePairInfo(groupIndex, -1, -1);
      ContactFacePairInfo& currentInfo = GetContactFacePairInfo(groupIndex, sideIndex, -1);

      if ( (currentInfo.typeIndex < generalInfo.typeIndex || generalInfo.typeIndex == IndexType(-1)) &&
         currentInfo.typeIndex != IndexType(-1) )
      {
        generalInfo.typeIndex = currentInfo.typeIndex;
        for(IndexType normalNumber = 0; normalNumber < 3; normalNumber++)
          generalInfo.normals[normalNumber] = Vector::zeroVector();
      }
      if(generalInfo.typeIndex == currentInfo.typeIndex)
      {
        for(IndexType normalNumber = 0; normalNumber < 3; normalNumber++)
          generalInfo.normals[normalNumber] += currentInfo.normals[normalNumber];
      }
    }
  }

  for(IndexType groupIndex = 0; groupIndex < GetContactFaceGroupsCount(); groupIndex++)
  {
    for(IndexType sideIndex0 = 0; sideIndex0 < GetContactFaceGroupSize(groupIndex) + 1; sideIndex0++)
    {
      for(IndexType sideIndex1 = 0; sideIndex1 < GetContactFaceGroupSize(groupIndex) + 1; sideIndex1++)
      {
        ContactFacePairInfo& currentInfo = GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1);
        if(currentInfo.typeIndex != IndexType(-1))
        {
          for(IndexType normalNumber = 0; normalNumber < 2; normalNumber++)
          {
            Scalar len = (currentInfo.normals[normalNumber]).Len();

            if(len != 0)
              currentInfo.normals[normalNumber] *= Scalar(1.0) / len;
            else
              currentInfo.normals[normalNumber] = Vector::zeroVector();
          }
        }
      }
    }
  }
}

Space3::IndexType GeomMesh<Space3>::
  GetContactFaceGroupsCount()
{
  return contactFaceGroupInfos.size();
}

Space3::IndexType GeomMesh<Space3>::
  GetContactFaceGroupSize(IndexType groupIndex)
{
  return contactFaceGroupInfos[groupIndex].size;
}

GeomMesh<Space3>::ContactFace GeomMesh<Space3>::
  GetContactFaceInGroup(IndexType groupIndex, IndexType faceIndex)
{
  return contactFaceGroupData[contactFaceGroupInfos[groupIndex].groupInfoOffset + faceIndex];
}

Space3::IndexType GeomMesh<Space3>::
  GetContactEdgeGroupsCount()
{
  return contactEdgeGroupInfos.size();
}

Space3::IndexType GeomMesh<Space3>::
  GetContactEdgeGroupSize(IndexType groupIndex)
{
  return contactEdgeGroupInfos[groupIndex].size;
}

GeomMesh<Space3>::ContactEdge GeomMesh<Space3>::
  GetContactEdgeInGroup(IndexType groupIndex, IndexType edgeIndex)
{
  return contactEdgeGroupData[contactEdgeGroupInfos[groupIndex].groupInfoOffset + edgeIndex];
}

Space3::IndexType GeomMesh<Space3>::
  GetContactNodeGroupsCount()
{
  return contactNodeGroupInfos.size();
}

Space3::IndexType GeomMesh<Space3>::
  GetContactNodeGroupSize(IndexType groupIndex)
{
  return contactNodeGroupInfos[groupIndex].size;
}

GeomMesh<Space3>::ContactNode GeomMesh<Space3>::
  GetContactNodeInGroup(IndexType groupIndex, IndexType nodeIndex)
{
  return contactNodeGroupData[contactNodeGroupInfos[groupIndex].groupInfoOffset + nodeIndex];
}

Space3::IndexType GeomMesh<Space3>::
  GetNodeSideInGroup(IndexType groupIndex, IndexType nodeIndex)
{
  for (IndexType sideIndex = 0; sideIndex < GetContactNodeGroupSize(groupIndex); sideIndex++)
  {
    if (GetContactNodeInGroup(groupIndex, sideIndex).nodeIndex == nodeIndex) return sideIndex;
  }
  return -1;
}

GeomMesh<Space3>::ContactNodePairInfo& GeomMesh<Space3>::
  GetContactNodePairInfo(IndexType groupIndex, IndexType sideIndex0, IndexType sideIndex1)
{
  if(sideIndex0 == IndexType(-1)) sideIndex0 = GetContactNodeGroupSize(groupIndex);
  if(sideIndex1 == IndexType(-1)) sideIndex1 = GetContactNodeGroupSize(groupIndex);

  return contactNodePairInfos[contactNodeGroupInfos[groupIndex].pairInfoOffset + sideIndex0 + 
    (contactNodeGroupInfos[groupIndex].size + 1) * sideIndex1];
}

Space3::IndexType GeomMesh<Space3>::
  GetEdgeSideInGroup(IndexType groupIndex, IndexType edgeIndex)
{
  for(IndexType sideIndex = 0; sideIndex < GetContactEdgeGroupSize(groupIndex); sideIndex++)
  {
    if(GetContactEdgeInGroup(groupIndex, sideIndex).edgeIndex == edgeIndex) return sideIndex;
  }
  return -1;
}

GeomMesh<Space3>::ContactEdgePairInfo& GeomMesh<Space3>::
  GetContactEdgePairInfo(IndexType groupIndex, IndexType sideIndex0, IndexType sideIndex1)
{
  if(sideIndex0 == IndexType(-1)) sideIndex0 = GetContactEdgeGroupSize(groupIndex);
  if(sideIndex1 == IndexType(-1)) sideIndex1 = GetContactEdgeGroupSize(groupIndex);

  return contactEdgePairInfos[contactEdgeGroupInfos[groupIndex].pairInfoOffset + sideIndex0 + 
    (contactEdgeGroupInfos[groupIndex].size + 1) * sideIndex1];
}

Space3::IndexType GeomMesh<Space3>::
  GetFaceSideInGroup(IndexType groupIndex, IndexType faceIndex)
{
  for(IndexType sideIndex = 0; sideIndex < GetContactFaceGroupSize(groupIndex); sideIndex++)
  {
    if(GetContactFaceInGroup(groupIndex, sideIndex).faceIndex == faceIndex) return sideIndex;
  }
  return -1;
}

GeomMesh<Space3>::ContactFacePairInfo& GeomMesh<Space3>::
  GetContactFacePairInfo(IndexType groupIndex, IndexType sideIndex0, IndexType sideIndex1)
{
  if(sideIndex0 == IndexType(-1)) sideIndex0 = GetContactFaceGroupSize(groupIndex);
  if(sideIndex1 == IndexType(-1)) sideIndex1 = GetContactFaceGroupSize(groupIndex);

  return contactFacePairInfos[contactFaceGroupInfos[groupIndex].pairInfoOffset + sideIndex0 + 
    (contactFaceGroupInfos[groupIndex].size + 1) * sideIndex1];
}

Space3::Vector GeomMesh<Space3>::GetFaceOuterNormal(IndexType *incidentNodes, Vector* vertexPositions)
{
  IndexType tempNodes[3];
  tempNodes[0] = incidentNodes[0];
  tempNodes[1] = incidentNodes[1];
  tempNodes[2] = incidentNodes[2];
  std::sort(tempNodes, tempNodes + 3);
  IndexType faceIndex = GetFaceIndex(tempNodes);

  if(faceIndex == IndexType(-1))
  {
    printf("wrong contact normal request\n");
    return Vector(-1, -1, -1);
  }

  return GetFaceOuterNormal(faceIndex, vertexPositions);
}

Space3::Vector GeomMesh<Space3>::
  GetFaceOuterNormal(IndexType faceIndex, Vector* vertexPositions)
{
  IndexType faceNode0 = faces[faceIndex].incidentNodes[0];
  IndexType faceNode1 = faces[faceIndex].incidentNodes[1];
  IndexType faceNode2 = faces[faceIndex].incidentNodes[2];

  IndexType cellIndex = faceTopologyInfos[faceIndex].incidentCells[0];

  IndexType cellNode0 = cells[cellIndex].incidentNodes[0];
  IndexType cellNode1 = cells[cellIndex].incidentNodes[1];
  IndexType cellNode2 = cells[cellIndex].incidentNodes[2];
  IndexType cellNode3 = cells[cellIndex].incidentNodes[3];

  IndexType internalPoint;
  if (cellNode0 != faceNode0 && cellNode0 != faceNode1 && cellNode0 != faceNode2)
  {
    internalPoint = cellNode0;
  } else
  if (cellNode1 != faceNode0 && cellNode1 != faceNode1 && cellNode1 != faceNode2)
  {
    internalPoint = cellNode1;
  } else
  if (cellNode2 != faceNode0 && cellNode2 != faceNode1 && cellNode2 != faceNode2)
  {
    internalPoint = cellNode2;
  } else
  if (cellNode3 != faceNode0 && cellNode3 != faceNode1 && cellNode3 != faceNode2)
  {
    internalPoint = cellNode3;
  } else
  {
    // fail
    internalPoint = cellNode0;
  }

  Vector origin;
  Vector vec0;
  Vector vec1;
  Vector vec2;

  if(vertexPositions != nullptr)
  {
    origin  = vertexPositions[faceNode0];
    vec0    = vertexPositions[faceNode1] - origin;
    vec1    = vertexPositions[faceNode2] - origin;
    vec2    = vertexPositions[internalPoint] - origin;
  } else
  {
    origin  = nodes[faceNode0].pos;
    vec0    = nodes[faceNode1].pos - origin;
    vec1    = nodes[faceNode2].pos - origin;
    vec2    = nodes[internalPoint].pos - origin;
  }

  if(MixedProduct(vec0, vec1, vec2) < 0.0)
  {
    return  vec0.GetNorm() ^ vec1.GetNorm(); //without normalization. weight is ~ proportional to sin(vec0 ^ vec1)
  }else
  {
    return -vec0.GetNorm() ^ vec1.GetNorm(); //without normalization. weight is ~ proportional to sin(vec0 ^ vec1)
  }
}

Space3::IndexType GeomMesh<Space3>::
  GetCellInNodeDirection(IndexType nodeIndex, Vector direction)
{
  Scalar surfaceEps = 0;//std::numeric_limits<Scalar>::epsilon();
  for(IndexType incidentCellIndex = 0; incidentCellIndex < nodeTopologyInfos[nodeIndex].incidentCellsCount; incidentCellIndex++)
  {
    IndexType cellIndex = nodeTopologyInfos[nodeIndex].incidentCells[incidentCellIndex];
    //CellNodeIterator nodeIterator = GetCellNodeIterator(cellIterator.GetCellIndex());

    Vector points[5]; //one last will be used just in case of overflow
    int pointsCount = 0;

    points[pointsCount++] = nodes[nodeIndex].pos;

    for(int i = 0; i < 4; i++)
    {
      if(nodeIndex != cells[cellIndex].incidentNodes[i]) //iterate all nodes != nodeIndex
      {
        points[pointsCount++] = nodes[cells[cellIndex].incidentNodes[i]].pos;
      }
    }

    assert(pointsCount == 4);

    for(int i = 1; i < 4; i++)
    {
      points[i] = points[i] - points[0];
    }

    Scalar side0 = MixedProduct(points[1], points[2], direction);
    Scalar side1 = MixedProduct(points[2], points[3], direction);
    Scalar side2 = MixedProduct(points[3], points[1], direction);
    Scalar side3 = MixedProduct(points[1], points[2], points[3]);

    if ((side0 >= -Scalar(surfaceEps)) && (side1 >= -Scalar(surfaceEps)) && (side2 >= -Scalar(surfaceEps)) && (side3 >= -Scalar(surfaceEps)))
    {
      return cellIndex;
    }

    if ((side0 <= Scalar(surfaceEps)) && (side1 <= Scalar(surfaceEps)) && (side2 <= Scalar(surfaceEps)) && (side3 <= Scalar(surfaceEps)))
    {
      return cellIndex;
    }
  };
  return -1;//there is no cell in this direction
};

Space3::IndexType GeomMesh<Space3>::
  GetCellInEdgeDirection(IndexType edgeIndex, Vector direction)
{
  Scalar surfaceEps = 0;//std::numeric_limits<Scalar>::epsilon();
  for(IndexType incidentCellIndex = 0; incidentCellIndex < edgeTopologyInfos[edgeIndex].incidentCellsCount; incidentCellIndex++)
  {
    IndexType cellIndex = edgeTopologyInfos[edgeIndex].incidentCells[incidentCellIndex];

    Vector points[5]; //one last will be used just in case of overflow
    IndexType pointsCount = 0;

    points[pointsCount++] = nodes[edges[edgeIndex].incidentNodes[0]].pos;
    points[pointsCount++] = nodes[edges[edgeIndex].incidentNodes[1]].pos;

    for(IndexType i = 0; i < 4; i++)
    {
       //iterate all nodes not incident to edges[edgeIndex]
      if((edges[edgeIndex].incidentNodes[0] != cells[cellIndex].incidentNodes[i]) &&
         (edges[edgeIndex].incidentNodes[1] != cells[cellIndex].incidentNodes[i]))
      {
        points[pointsCount++] = nodes[cells[cellIndex].incidentNodes[i]].pos;
      }
    }

    assert(pointsCount == 4);

    for(IndexType i = 1; i < 4; i++)
    {
      points[i] = points[i] - points[0];
    }

    Scalar side0 = MixedProduct(points[1], points[2], direction);
    Scalar side1 = MixedProduct(points[3], points[1], direction);
    Scalar side2 = MixedProduct(points[1], points[2], points[3]);

    if ((side0 >= -Scalar(surfaceEps)) && (side1 >= -Scalar(surfaceEps)) && (side2 >= -Scalar(surfaceEps)))
    {
      return cellIndex;
    }

    if ((side0 <= Scalar(surfaceEps)) && (side1 <= Scalar(surfaceEps)) && (side2 <= Scalar(surfaceEps)))
    {
      return cellIndex;
    }

    /*if(((side0 * side1) >= 0.0) && ((side1 * side2) >= 0.0)) //sgn(side0) == sgn(side1) == sgn(side2)
    {
      return cellIndex;
    }*/
  };
  return -1;//there is no cell in this direction
}

Space3::IndexType GeomMesh<Space3>::
  GetCellInFaceDirection(IndexType faceIndex, Vector direction)
{
  Scalar surfaceEps = 0;//std::numeric_limits<Scalar>::epsilon();
  for(IndexType incidentCellIndex = 0; incidentCellIndex < faceTopologyInfos[faceIndex].incidentCellsCount; incidentCellIndex++)
  {
    IndexType cellIndex = faceTopologyInfos[faceIndex].incidentCells[incidentCellIndex];

    Vector points[3];
    for(int i = 0; i < 3; i++)
    {
      points[i] = nodes[faces[faceIndex].incidentNodes[i]].pos;
    }

    for(int i = 1; i < 3; i++)
    {
      points[i] = points[i] - points[0]; //move reference frame to points[0]
    }

    Scalar referenceSide = MixedProduct(points[1], points[2], direction);

    for(IndexType i = 0; i < 4; i++)
    {
      if((faces[faceIndex].incidentNodes[0] != cells[cellIndex].incidentNodes[i]) &&
         (faces[faceIndex].incidentNodes[1] != cells[cellIndex].incidentNodes[i]) &&
         (faces[faceIndex].incidentNodes[2] != cells[cellIndex].incidentNodes[i]))
      {
        Vector testPoint = nodes[cells[cellIndex].incidentNodes[i]].pos - points[0];

        Scalar side = MixedProduct(points[1], points[2], testPoint);
        if(referenceSide * side >= -Scalar(surfaceEps))
        {
          return cellIndex;
        }
      }
    }
  };
  return -1;//there is no cell in this direction
}

Space3::IndexType GeomMesh<Space3>::
  GetCellIndex(IndexType node0, IndexType node1, IndexType node2, IndexType node3) const
{
  IndexType nodeIndices[Space::NodesPerCell];
  nodeIndices[0] = node0;
  nodeIndices[1] = node1;
  nodeIndices[2] = node2;
  nodeIndices[3] = node3;
  return GetCellIndex(nodeIndices);
}

Space3::IndexType GeomMesh<Space3>::
  GetFaceIndex(const IndexType nodeIndices[Space::NodesPerFace]) const
{
  // node indices have to be ordered
  IndexType orderedIndices[Space::NodesPerFace];
  std::copy(nodeIndices, nodeIndices + Space::NodesPerFace, orderedIndices);
  std::sort(orderedIndices, orderedIndices + Space::NodesPerFace);

  IndexType baseNode = orderedIndices[0]; //may be any of given nodes

  for(IndexType incidentCellIndex = 0; incidentCellIndex < nodeTopologyInfos[baseNode].incidentCellsCount; incidentCellIndex++)
  {
    IndexType cellIndex = nodeTopologyInfos[baseNode].incidentCells[incidentCellIndex];
    for(IndexType incidentFaceIndex = 0; incidentFaceIndex < Space::FacesPerCell; incidentFaceIndex++)
    {
      IndexType faceIndex = cellTopologyInfos[cellIndex].incidentFaces[incidentFaceIndex];

      bool found = true;
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
      {
        if (faces[faceIndex].incidentNodes[nodeNumber] != orderedIndices[nodeNumber])
        {
          found = false;
          break;
        }
      }
      if (found)
      {
        return faceIndex;
      }
    }
  }
  return IndexType(-1); //there is no such face
}

Space3::IndexType GeomMesh<Space3>::
  GetFaceIndex(IndexType node0, IndexType node1, IndexType node2) const
{
  IndexType nodeIndices[Space::NodesPerFace];
  nodeIndices[0] = node0;
  nodeIndices[1] = node1;
  nodeIndices[2] = node2;
  return GetFaceIndex(nodeIndices);
}

Space3::IndexType GeomMesh<Space3>::GetSubmeshNodesCount(IndexType submeshIndex)
{
  return submeshInfos[submeshIndex].nodesCount;
}

Space3::IndexType GeomMesh<Space3>::GetSubmeshEdgesCount(IndexType submeshIndex)
{
  return submeshInfos[submeshIndex].edgesCount;
}

GeomMesh<Space3>::Face* GeomMesh<Space3>::GetFace(IndexType index)
{
  return &(faces[index]);
}

const GeomMesh<Space3>::Face* GeomMesh<Space3>::GetFace(IndexType index) const
{
  return &(faces[index]);
}

Space3::IndexType GeomMesh<Space3>::GetFacesCount()
{
  return IndexType(faces.size());
}

Space3::IndexType GeomMesh<Space3>::GetSubmeshFacesCount(IndexType submeshIndex)
{
  return submeshInfos[submeshIndex].facesCount;
}

Space3::IndexType GeomMesh<Space3>::GetSubmeshCellsCount(IndexType submeshIndex)
{
  return IndexType(submeshInfos[submeshIndex].cellsCount);
}

Space3::IndexType GeomMesh<Space3>::GetSubmeshesCount()
{
  return submeshesCount;
}

Space3::IndexType GetRelativeOrientationSide(Space3::IndexType* firstNodes, Space3::IndexType* secondNodes)
{
  Space3::IndexType i, j, k;
  Space3::IndexType a, b, c;
  Space3::IndexType R;

  if (firstNodes[0]<firstNodes[1])
  {
    if (firstNodes[1]<firstNodes[2])
    {
      //012
      i=0;
      j=1;
      k=2;
    }
    else
    {
      if (firstNodes[0]<firstNodes[2])
      {
        //021
        i=0;
        j=2;
        k=1;
      }
      else
      {
        //120
        i=2;
        j=0;
        k=1;
      };
    };
  }
  else
  {
    if (firstNodes[1]>firstNodes[2])
    {
      //210
      i=2;
      j=1;
      k=0;
    }
    else
    {
      if (firstNodes[0]<firstNodes[2])
      {
        //102
        i=1;
        j=0;
        k=2;
      }
      else
      {
        //201
        i=1;
        j=2;
        k=0;
      };
    };
  };  

  a=secondNodes[i];
  b=secondNodes[j];
  c=secondNodes[k];

  if (a<b)
  {
    if (b<c)
    {
      R=12;
    }
    else
    {
      if (a<c)
      {
        R=1021;
      }
      else
      {
        R=120;
      };
    };
  }
  else
  {
    if (b>c)
    {
      R=1210;
    }
    else
    {
      if (a<c)
      {
        R=1102;
      }
      else
      {
        R=201;
      };
    };
  };  

  return(R);
}

Space3::IndexType GetRelativeOrientationEdge(Space3::IndexType* firstNodes, Space3::IndexType* secondNodes)
{
  int R;
  if (firstNodes[0]<firstNodes[1])
  {
    if (secondNodes[0]<secondNodes[1])
    {
      R=0;
    }
    else
    {
      R=1;
    };
  }
  else
  {
    if (secondNodes[0]<secondNodes[1])
    {
      R=1;
    }
    else
    {
      R=0;
    };
  };
  return(R);
}

int GetCorrespondingFaceIndex(int relativeOrientation, int firstIndex, int number)
{
  int L;
  int P[3];
  int N;
  int R;

  N=relativeOrientation;
  L=N/1000;

  N=N-1000*L;
  P[0]=N/100;

  N=N-100*P[0];
  P[1]=N/10;

  N=N-10*P[1];
  P[2]=N;

  if (number<4)
  {
    R=firstIndex;
  }

  else if (number==4)
  {
    R=P[firstIndex];
  }

  else if (number==5)
  {
    if (firstIndex<3)
    {
      R=P[firstIndex];
    }
    else
    {
      if (L==0)
      {
        R=P[firstIndex-3]+3;
      }
      else
      {
        if (firstIndex==3)
        {
          R=P[1]+3;
        }
        else if (firstIndex==4)
        {
          R=P[2]+3;
        }
        else if (firstIndex==5)
        {
          R=P[0]+3;
        };
      };
    };
  };
  return(R);
}

int GetCorrespondingEdgeIndex(int relativeOrientation, int firstIndex, int number)
{
  int R;
  if (number > 2 && relativeOrientation == 1)
  {
    R = number - 2 - firstIndex;
  } else 
  {
    R = firstIndex;
  }
  return R;
}

bool PointIsInCell(Space3::Vector cellPoint0, Space3::Vector cellPoint1, 
  Space3::Vector cellPoint2, Space3::Vector cellPoint3, Space3::Vector testPoint)
{
  typedef Space3::Scalar Scalar;
  //012
  Scalar side0 = MixedProduct(cellPoint1 - cellPoint0, cellPoint2 - cellPoint0, cellPoint3 - cellPoint0) *
                 MixedProduct(cellPoint1 - cellPoint0, cellPoint2 - cellPoint0, testPoint  - cellPoint0);
  //013
  Scalar side1 = MixedProduct(cellPoint1 - cellPoint0, cellPoint3 - cellPoint0, cellPoint2 - cellPoint0) *
                 MixedProduct(cellPoint1 - cellPoint0, cellPoint3 - cellPoint0, testPoint  - cellPoint0);
  //023
  Scalar side2 = MixedProduct(cellPoint2 - cellPoint0, cellPoint3 - cellPoint0, cellPoint1 - cellPoint0) *
                 MixedProduct(cellPoint2 - cellPoint0, cellPoint3 - cellPoint0, testPoint  - cellPoint0);
  //123
  Scalar side3 = MixedProduct(cellPoint2 - cellPoint1, cellPoint3 - cellPoint1, cellPoint0 - cellPoint1) *
                 MixedProduct(cellPoint2 - cellPoint1, cellPoint3 - cellPoint1, testPoint  - cellPoint1);

  return ((side0 >= 0.0) && (side1 >= 0.0) && (side2 >= 0.0) && (side3 >= 0.0));
}

void GeomMesh<Space3>::BuildCellAdditionalTopology(IndexType *internalContactTypes)
{
  additionalCellInfos.resize(GetCellsCount());
  for(IndexType cellIndex = 0; cellIndex < GetCellsCount(); cellIndex++)
  {
    for(IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
    {
      additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex = IndexType(-1);
      additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber = IndexType(-1);
      additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType = IndexType(-1);
      additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].orientation = IndexType(-1);
    }
  }

  for(IndexType faceIndex = 0; faceIndex < GetFacesCount(); faceIndex++)
  {
    IndexType srcFaceNumber;
    IndexType dstFaceNumber;
    IndexType orientation;

    if(faceTopologyInfos[faceIndex].incidentCellsCount == 2)
    {
      IndexType contactingCells[2];
      for(IndexType contactingCellIndex = 0; contactingCellIndex < 2; contactingCellIndex++)
      {
        contactingCells[contactingCellIndex] = faceTopologyInfos[faceIndex].incidentCells[contactingCellIndex];
      }

      for(IndexType cellNumber = 0; cellNumber < 2; cellNumber++)
      {
        IndexType srcCellIndex = contactingCells[cellNumber];
        IndexType dstCellIndex = contactingCells[(cellNumber + 1) % 2];

        if(FindCommonFace(srcCellIndex, dstCellIndex, srcFaceNumber, dstFaceNumber, orientation))
        {
          additionalCellInfos[srcCellIndex].neighbouringFaces[srcFaceNumber].correspondingCellIndex = dstCellIndex;
          additionalCellInfos[srcCellIndex].neighbouringFaces[srcFaceNumber].correspondingFaceNumber = dstFaceNumber;
          if(internalContactTypes)
            additionalCellInfos[srcCellIndex].neighbouringFaces[srcFaceNumber].interactionType = std::max(internalContactTypes[srcCellIndex], internalContactTypes[dstCellIndex]); //priority is given to a larger type;
          else
            additionalCellInfos[srcCellIndex].neighbouringFaces[srcFaceNumber].interactionType = IndexType(0);
          additionalCellInfos[srcCellIndex].neighbouringFaces[srcFaceNumber].orientation = orientation;
        } else
        {
          assert(0); //topology error
        }
      }
    }
  }
}

Space3::Scalar GeomMesh<Space3>::GetFaceSquare(Vector faceVertices[Space::NodesPerFace])
{
  return Scalar(0.5) * ((faceVertices[2] - faceVertices[0]) ^ (faceVertices[1] - faceVertices[0])).Len();
}

Space3::IndexType GeomMesh<Space3>::GetFacePairOrientation(IndexType srcCellNodes[3], IndexType dstCellNodes[3])
{
  if(srcCellNodes[0] == dstCellNodes[0] && srcCellNodes[1] == dstCellNodes[2] && srcCellNodes[2] == dstCellNodes[1])
    return 0;
  if(srcCellNodes[0] == dstCellNodes[1] && srcCellNodes[1] == dstCellNodes[0] && srcCellNodes[2] == dstCellNodes[2])
    return 1;
  if(srcCellNodes[0] == dstCellNodes[2] && srcCellNodes[1] == dstCellNodes[1] && srcCellNodes[2] == dstCellNodes[0])
    return 2;
  return IndexType(-1);
}

bool GeomMesh<Space3>::FindCommonFace(IndexType srcCellIndex, IndexType dstCellIndex, 
  IndexType &_srcFaceNumber, IndexType &_dstFaceNumber, IndexType &_orientation)
{
  IndexType srcCellNodes[4];
  IndexType dstCellNodes[4];

  GetFixedCellIndices(srcCellIndex, srcCellNodes);
  GetFixedCellIndices(dstCellIndex, dstCellNodes);

  for(IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
  {
    IndexType srcFaceNodes[3];
    GetCellFaceNodes(srcCellNodes, srcFaceNumber, srcFaceNodes);

    for(IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
    {
      IndexType dstFaceNodes[3];
      GetCellFaceNodes(dstCellNodes, dstFaceNumber, dstFaceNodes);

      IndexType orientation = GetFacePairOrientation(srcFaceNodes, dstFaceNodes);

      if(orientation != IndexType(-1))
      {
        _orientation = orientation;
        _srcFaceNumber = srcFaceNumber;
        _dstFaceNumber = dstFaceNumber;
        return 1;
      }
    }
  }
  return 0;
}

Space3::IndexType
GeomMesh<Space3>::GetFaceNumber(IndexType cellIndex, const IndexType faceNodeIndices[Space::NodesPerFace]) const
{
  IndexType orderedIndices[Space::NodesPerFace];
  std::copy(faceNodeIndices, faceNodeIndices + Space::NodesPerFace, orderedIndices);
  std::sort(orderedIndices, orderedIndices + Space::NodesPerFace);

  IndexType cellIncidentNodes[Space::NodesPerCell];
  GetFixedCellIndices(cellIndex, cellIncidentNodes);
  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
  {
    IndexType faceIncidentNodes[Space::NodesPerFace];
    GetCellFaceNodes(cellIncidentNodes, faceNumber, faceIncidentNodes);
    std::sort(faceIncidentNodes, faceIncidentNodes + Space::NodesPerFace);
    bool found = true;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
    {
      if (orderedIndices[nodeNumber] != faceIncidentNodes[nodeNumber]) found = false;
    }
    if (found) return faceNumber;
  }
  return IndexType(-1);
}

Space3::IndexType 
GeomMesh<Space3>::GetFaceNumber(IndexType cellIndex, IndexType nodeIndex0, IndexType nodeIndex1, IndexType nodeIndex2) const
{
  IndexType faceNodeIndices[Space::NodesPerFace];
  faceNodeIndices[0] = nodeIndex0;
  faceNodeIndices[1] = nodeIndex1;
  faceNodeIndices[2] = nodeIndex2;
  return GetFaceNumber(cellIndex, faceNodeIndices);
}

GeomMesh<Space3>::FaceLocationPair GeomMesh<Space3>::GetFaceLocation(IndexType nodeIndex0, IndexType nodeIndex1, IndexType nodeIndex2) const
{
  IndexType nodeIndices[Space::NodesPerFace];
  nodeIndices[0] = nodeIndex0;
  nodeIndices[1] = nodeIndex1;
  nodeIndices[2] = nodeIndex2;
  return GetFaceLocation(nodeIndices);
}

GeomMesh<Space3>::FaceLocationPair GeomMesh<Space3>::GetFaceLocation(const IndexType nodeIndices[Space::NodesPerFace]) const
{
  FaceLocationPair faceLocationPair;

  IndexType faceIndex = GetFaceIndex(nodeIndices);
  assert(faceIndex != IndexType(-1));

  for (IndexType incidentCellNumber = 0; incidentCellNumber < faceTopologyInfos[faceIndex].incidentCellsCount; ++incidentCellNumber)
  {
    IndexType cellIndex = faceTopologyInfos[faceIndex].incidentCells[incidentCellNumber];
    faceLocationPair.faces[incidentCellNumber].cellIndex = cellIndex;
    IndexType faceNumber = GetFaceNumber(cellIndex, nodeIndices);
    assert(faceNumber != IndexType(-1));
    faceLocationPair.faces[incidentCellNumber].faceNumber = faceNumber;
  }

  return faceLocationPair;
}

void GeomMesh<Space3>::GetCellFaceNodes(IndexType cellIndex, IndexType faceNumber, IndexType* faceNodes) const
{
  IndexType incidentNodes[Space::NodesPerCell];
  GetFixedCellIndices(cellIndex, incidentNodes);
  GetCellFaceNodes(incidentNodes, faceNumber, faceNodes);
}

void GeomMesh<Space3>::GetCellFaceNodes(const IndexType *cellIncidentNodes, IndexType faceNumber, IndexType *faceNodes) const
{
  switch(faceNumber)
  {
    case 0:
    {
      faceNodes[0] = cellIncidentNodes[0];
      faceNodes[1] = cellIncidentNodes[2];
      faceNodes[2] = cellIncidentNodes[1];
    } break;
    case 1:
    {
      faceNodes[0] = cellIncidentNodes[0];
      faceNodes[1] = cellIncidentNodes[1];
      faceNodes[2] = cellIncidentNodes[3];
    } break;
    case 2:
    {
      faceNodes[0] = cellIncidentNodes[0];
      faceNodes[1] = cellIncidentNodes[3];
      faceNodes[2] = cellIncidentNodes[2];
    } break;
    case 3:
    {
      faceNodes[0] = cellIncidentNodes[1];
      faceNodes[1] = cellIncidentNodes[2];
      faceNodes[2] = cellIncidentNodes[3];
    } break;
  }
}

void GeomMesh<Space3>::
  BuildBoundariesInfo(BoundaryFace* boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount)
{
  IndexType boundaryFaceIndex = 0;
  for(IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; boundaryTypeIndex++)
  {
    for(IndexType typeFaceIndex = 0; typeFaceIndex < boundaryFacesCount[boundaryTypeIndex]; typeFaceIndex++)
    {
      IndexType faceIndex;
      IndexType cellIndex;
      IndexType localFaceNumber;

      IndexType incidentNodes[3];
      for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
      {
        incidentNodes[nodeNumber] = boundaryFaces[boundaryFaceIndex].nodeIndices[nodeNumber];
      }
      faceIndex = GetFaceIndex(incidentNodes);
      if (faceIndex == IndexType(-1)) continue;

      for(IndexType sideNumber = 0; sideNumber < faceTopologyInfos[faceIndex].incidentCellsCount; sideNumber++)
      {
        cellIndex = faceTopologyInfos[faceIndex].incidentCells[sideNumber];
        assert(cellIndex != IndexType(-1));

        localFaceNumber = 
          GetLocalFaceNumber(cellIndex, boundaryFaces[boundaryFaceIndex].nodeIndices);

        assert(localFaceNumber != IndexType(-1));

        AdditionalCellInfo<Space3>::NeighbourFace& info =
          additionalCellInfos[cellIndex].neighbouringFaces[localFaceNumber];

        info.orientation = 0;
        info.correspondingCellIndex   = IndexType(-1);
        info.correspondingFaceNumber  = IndexType(-1);
        info.interactionType          = boundaryTypeIndex;
      }
      boundaryFaceIndex++;
    }
  }
}

Space3::IndexType GeomMesh<Space3>::
  GetLocalFaceNumber(IndexType cellIndex, IndexType* faceIndices)
{
  IndexType cellNodeIndices[4];
  GetFixedCellIndices(cellIndex, cellNodeIndices);

  for(IndexType localFaceNumber = 0; localFaceNumber < 4; localFaceNumber++)
  {
    IndexType testFaceNodes[3];
    GetCellFaceNodes(cellNodeIndices, localFaceNumber, testFaceNodes);

    IndexType matches = 0;
    for(IndexType nodeNumber0 = 0; nodeNumber0 < 3; nodeNumber0++)
    {
      for(IndexType nodeNumber1 = 0; nodeNumber1 < 3; nodeNumber1++)
      {
        if(testFaceNodes[nodeNumber0] == faceIndices[nodeNumber1]) matches++;
      }
    }
    if(matches == 3) return localFaceNumber;
  }
  return IndexType(-1);
}

void GeomMesh<Space3>::BuildContactsInfo(FacePairIndices* contactFaces, 
  IndexType* contactFacesCount, IndexType contactTypesCount)
{
  IndexType rejectedContacts = 0;

  IndexType contactFaceIndex = 0;
  for(IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; contactTypeIndex++)
  {
    for(IndexType typeFaceIndex = 0; typeFaceIndex < contactFacesCount[contactTypeIndex]; typeFaceIndex++)
    {
      IndexType faceIndices[2];
      IndexType cellIndices[2];
      IndexType localFaceNumbers[2];

      bool fineContact = true;

      for(IndexType sideNumber = 0; sideNumber < 2; sideNumber++)
      {
        IndexType incidentNodes[Space::NodesPerFace];
        for(IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; nodeNumber++)
        {
          incidentNodes[nodeNumber] = contactFaces[contactFaceIndex].faces[sideNumber].nodeIndices[nodeNumber];
        }
        faceIndices[sideNumber] = GetFaceIndex(incidentNodes);
        if(faceIndices[sideNumber] == IndexType(-1))
        {
          fineContact = false;
          continue;
        }
      }

      if(!fineContact)
      {
        rejectedContacts++;
        continue;
      }

      if(faceTopologyInfos[faceIndices[0]].incidentCellsCount == 1 && faceTopologyInfos[faceIndices[1]].incidentCellsCount == 1)
      {
        for(IndexType sideNumber = 0; sideNumber < 2; sideNumber++)
        {
          assert(faceTopologyInfos[faceIndices[sideNumber]].incidentCellsCount == 1);

          cellIndices[sideNumber] = faceTopologyInfos[faceIndices[sideNumber]].incidentCells[0];
          assert(cellIndices[sideNumber] != -1);

          localFaceNumbers[sideNumber] = 
            GetLocalFaceNumber(cellIndices[sideNumber], contactFaces[contactFaceIndex].faces[sideNumber].nodeIndices);
          assert(localFaceNumbers[sideNumber] != -1);
        }
      }else
      if(faceIndices[0] == faceIndices[1] && faceTopologyInfos[faceIndices[0]].incidentCellsCount == 2)
      {
        IndexType commonFace = faceIndices[0];
        for(IndexType sideNumber = 0; sideNumber < 2; sideNumber++)
        {
          cellIndices[sideNumber] = faceTopologyInfos[commonFace].incidentCells[sideNumber];
          assert(cellIndices[sideNumber] != -1);

          localFaceNumbers[sideNumber] = 
            GetLocalFaceNumber(cellIndices[sideNumber], contactFaces[contactFaceIndex].faces[sideNumber].nodeIndices);
          assert(localFaceNumbers[sideNumber] != -1);
        }
      }

      for(IndexType sideNumber = 0; sideNumber < 2; sideNumber++)
      {
        IndexType srcSide = sideNumber;
        IndexType dstSide = (sideNumber + 1) % 2;

        IndexType srcCellNodes[Space::NodesPerCell];
        IndexType dstCellNodes[Space::NodesPerCell];

        GetFixedCellIndices(cellIndices[srcSide], srcCellNodes);
        GetFixedCellIndices(cellIndices[dstSide], dstCellNodes);

        IndexType srcFaceNodes[3];
        IndexType dstFaceNodes[3];

        GetCellFaceNodes(srcCellNodes, localFaceNumbers[srcSide], srcFaceNodes);
        GetCellFaceNodes(dstCellNodes, localFaceNumbers[dstSide], dstFaceNodes);

        IndexType srcFaceRelativeNodes[3];
        IndexType dstFaceRelativeNodes[3];


        for(IndexType contactNodeNumber = 0; contactNodeNumber < 3; contactNodeNumber++)
        {
          IndexType srcNodeNumber = IndexType(-1);
          IndexType dstNodeNumber = IndexType(-1);

          for(IndexType localNodeNumber = 0; localNodeNumber < 3; localNodeNumber++)
          {
            if(contactFaces[contactFaceIndex].faces[srcSide].nodeIndices[contactNodeNumber] == srcFaceNodes[localNodeNumber])
            {
              srcNodeNumber = localNodeNumber;
            }
          }

          for(IndexType localNodeNumber = 0; localNodeNumber < 3; localNodeNumber++)
          {
            if(contactFaces[contactFaceIndex].faces[dstSide].nodeIndices[contactNodeNumber] == dstFaceNodes[localNodeNumber])
            {
              dstNodeNumber = localNodeNumber;
            }
          }
          assert(srcNodeNumber != -1);
          assert(dstNodeNumber != -1);
          dstFaceRelativeNodes[srcNodeNumber] = dstNodeNumber;
        }
        for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
        {
          srcFaceRelativeNodes[nodeNumber] = nodeNumber;
          assert(dstFaceRelativeNodes[nodeNumber] != IndexType(-1));
        }

        AdditionalCellInfo<Space3>::NeighbourFace &info =
          additionalCellInfos[cellIndices[srcSide]].neighbouringFaces[localFaceNumbers[srcSide]];

        info.orientation              = GetFacePairOrientation(srcFaceRelativeNodes, dstFaceRelativeNodes);
        assert(info.orientation != IndexType(-1));
        if(info.orientation != IndexType(-1))
        {
          info.correspondingCellIndex   = cellIndices     [dstSide];
          info.correspondingFaceNumber  = localFaceNumbers[dstSide];
          info.interactionType          = contactTypeIndex;
        }
      }
      contactFaceIndex++;
    }
  }

  if(rejectedContacts != 0)
  {
    printf("Volume mesh has rejected %d contact faces\n", (int)rejectedContacts);
  }
}

void GeomMesh<Space3>::BuildAdditionalTopology(
  FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
  BoundaryFace*    boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount,
  IndexType *internalContactTypes)
{
  printf("Building geom mesh additional topology \n");
  BuildCellAdditionalTopology(internalContactTypes);
  BuildContactsInfo(contactFaces, contactFacesCount, contactTypesCount);
  BuildBoundariesInfo(boundaryFaces, boundaryFacesCount, boundaryTypesCount);

  if (internalContactTypes)
  {
    for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
    {
      for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
      {
        IndexType correspondingCellIndex = GetCorrespondingCellIndex(cellIndex, faceNumber);
        if (correspondingCellIndex != IndexType(-1))
        {
          IndexType interactionType = std::max(internalContactTypes[cellIndex], internalContactTypes[correspondingCellIndex]);
          additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType = std::max(interactionType,
            additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType);
        }
      }
    }
  }
}

void GeomMesh<Space3>::GetCellFaceVertices(IndexType cellIndex, IndexType faceNumber, Vector* faceVertices) const
{
  IndexType faceNodeIndices[Space::NodesPerFace];
  GetCellFaceNodes(cellIndex, faceNumber, faceNodeIndices);
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; nodeNumber++)
  {
    faceVertices[nodeNumber] = nodes[faceNodeIndices[nodeNumber]].pos;
  }
}

Space3::Vector GeomMesh<Space3>::GetFaceExternalNormal(Vector* faceGlobalVertices) const
{
  Vector externalNormal = (faceGlobalVertices[1] - faceGlobalVertices[0]) ^ (faceGlobalVertices[2] - faceGlobalVertices[0]);
  return externalNormal.GetNorm();
}

Space3::Vector GeomMesh<Space3>::GetFaceExternalNormal(IndexType cellIndex, IndexType faceNumber) const
{
  Vector faceGlobalVertices[Space::NodesPerFace];
  GetCellFaceVertices(cellIndex, faceNumber, faceGlobalVertices);
  return GetFaceExternalNormal(faceGlobalVertices);
}

void GeomMesh<Space3>::GetGhostCellVertices(IndexType cellIndex, IndexType boundaryFaceNumber, Vector* ghostCellVertices) const
{
  Vector points[Space::NodesPerCell];
  GetCellVertices(cellIndex, points);

  for (IndexType vertexNumber = 0; vertexNumber < Space::NodesPerCell; ++vertexNumber)
  {
    ghostCellVertices[vertexNumber] = points[vertexNumber];
  }

  IndexType mirroringNodeNumber = 3 - boundaryFaceNumber;
  Vector faceNormal = GetFaceExternalNormal(cellIndex, boundaryFaceNumber);

  ghostCellVertices[mirroringNodeNumber] = points[mirroringNodeNumber] -
    faceNormal * (faceNormal * (points[mirroringNodeNumber] - points[boundaryFaceNumber])) * Scalar(2.0);

  std::swap(ghostCellVertices[0], ghostCellVertices[1]);
}

GeomMesh<Space3>::FaceLocationPair GeomMesh<Space3>::GetFaceLocation(const FacePairIndices& facePairIndices) const
{  
  FaceLocationPair faceLocationPair;
  IndexType facesFoundCount = 0;

  for (IndexType facePairIndex = 0; facePairIndex < 2; ++facePairIndex)
  {
    FaceLocationPair facesLocation = GetFaceLocation(facePairIndices.faces[facePairIndex].nodeIndices);

    for (IndexType faceIndex = 0; faceIndex < 2; ++faceIndex)
    {
      if (!facesLocation.faces[faceIndex].IsNull())
      {
        bool found = false;
        for (IndexType facesFoundIndex = 0; facesFoundIndex < facesFoundCount; ++facesFoundIndex)
        {
          if (faceLocationPair.faces[facesFoundIndex] == facesLocation.faces[faceIndex]) found = true;
        }
        if (!found) faceLocationPair.faces[facesFoundCount++] = facesLocation.faces[faceIndex];
      }
    }
  }
  assert(facesFoundCount == 2);
  return faceLocationPair;
}

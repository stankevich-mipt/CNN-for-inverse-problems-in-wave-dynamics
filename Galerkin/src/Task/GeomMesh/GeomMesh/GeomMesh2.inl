GeomMesh<Space2>::BoundaryEdge::BoundaryEdge(): EdgeIndices()
{
}

GeomMesh<Space2>::BoundaryEdge::BoundaryEdge(IndexType nodeIndex0, IndexType nodeIndex1): 
  EdgeIndices(nodeIndex0, nodeIndex1)
{
}

void GeomMesh<Space2>::BuildAdditionalTopology(
  EdgePairIndices* contactEdges,  IndexType* contactEdgesCount,  IndexType contactTypesCount,
  BoundaryEdge*    boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount,
  IndexType *internalContactTypes)
{
  printf("Building geom mesh additional topology \n");
  BuildCellAdditionalTopology(internalContactTypes);
  BuildContactsInfo(contactEdges, contactEdgesCount, contactTypesCount);
  BuildBoundariesInfo(boundaryEdges, boundaryEdgesCount, boundaryTypesCount);

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
          additionalCellInfos[cellIndex].neighbouringEdges[faceNumber].interactionType = std::max(interactionType,
            additionalCellInfos[cellIndex].neighbouringEdges[faceNumber].interactionType);
        }
      }
    }
  }
}

GeomMesh<Space2>::EdgeLocationPair
GeomMesh<Space2>::GetEdgeLocation(IndexType nodeIndex0, IndexType nodeIndex1) const
{
  EdgeLocationPair edgeLocationPair;

  IndexType edgeLocationCount = 0;

  for (IndexType incidentCellNumber = 0; incidentCellNumber < nodeTopologyInfos[nodeIndex0].incidentCellsCount; ++incidentCellNumber) 
  {
    IndexType incidentCellIndex = nodeTopologyInfos[nodeIndex0].incidentCells[incidentCellNumber];
    IndexType edgeNumber = GetEdgeNumber(incidentCellIndex , nodeIndex0, nodeIndex1);
    if (edgeNumber != IndexType(-1)) 
    {
      edgeLocationPair.edges[edgeLocationCount] = EdgeLocation(incidentCellIndex, edgeNumber);
      ++edgeLocationCount;
    }
  }

  return edgeLocationPair;
}

GeomMesh<Space2>::EdgeLocationPair
GeomMesh<Space2>::GetEdgeLocation(EdgeIndices edgeIndices) const
{
  return GetEdgeLocation(edgeIndices.nodeIndices[0], edgeIndices.nodeIndices[1]);
}

Space2::IndexType 
GeomMesh<Space2>::GetEdgeNumber(IndexType cellIndex, IndexType nodeIndex0, IndexType nodeIndex1) const
{
  IndexType cellIncidentNodes[3];
  GetFixedCellIndices(cellIndex, cellIncidentNodes);
  for (IndexType nodeNumber = 0; nodeNumber < 3; ++nodeNumber) 
  {
    IndexType edgeIncidentNodes[2];
    GetCellEdgeNodes(cellIncidentNodes, nodeNumber, edgeIncidentNodes);
  
    if ( (edgeIncidentNodes[0] == nodeIndex0 && edgeIncidentNodes[1] == nodeIndex1) || 
         (edgeIncidentNodes[0] == nodeIndex1 && edgeIncidentNodes[1] == nodeIndex0) )
    {
      return nodeNumber;
    }
  }
  return IndexType(-1);
}

GeomMesh<Space2>::EdgeLocationPair
GeomMesh<Space2>::BuildEdgeLocation(IndexType nodeIndex1, IndexType nodeIndex2)
{
  EdgeLocationPair edgeLocationPair;
  if (nodeIndex1 != IndexType(-1) && nodeIndex2 != IndexType(-1)) 
  {
    edgeLocationPair = GetEdgeLocation(nodeIndex1, nodeIndex2);
  }
  return edgeLocationPair;
}

void GeomMesh<Space2>::BuildContactsInfo(EdgePairIndices* contactEdges, IndexType* contactEdgesCount, IndexType contactTypesCount)
{
  IndexType offset = 0;

  for (IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; ++contactTypeIndex)
  {    
    for (IndexType contactEdgeIndex = 0; contactEdgeIndex < contactEdgesCount[contactTypeIndex]; ++contactEdgeIndex)
    {
      EdgePairIndices edgePairIndices = contactEdges[offset + contactEdgeIndex];
      
      EdgeLocationPair contactEdgesLocation[2];

      for (IndexType edgePairIndex = 0; edgePairIndex < 2; ++edgePairIndex)
      {
        IndexType nodeIndex1 = edgePairIndices.edges[edgePairIndex].nodeIndices[0];
        IndexType nodeIndex2 = edgePairIndices.edges[edgePairIndex].nodeIndices[1];
        contactEdgesLocation[edgePairIndex] = BuildEdgeLocation(nodeIndex1, nodeIndex2);
      }
  
      for (IndexType edgeIndex0 = 0; edgeIndex0 < 2; ++edgeIndex0)
      {
        for (IndexType edgeIndex1 = 0; edgeIndex1 < 2; ++edgeIndex1)
        {                  
          if (!contactEdgesLocation[0].edges[edgeIndex0].IsNull() && !contactEdgesLocation[1].edges[edgeIndex1].IsNull() &&
              contactEdgesLocation[0].edges[edgeIndex0] != contactEdgesLocation[1].edges[edgeIndex1])
          {
            IndexType cellIndex0 = contactEdgesLocation[0].edges[edgeIndex0].cellIndex;
            IndexType edgeNumber0 = contactEdgesLocation[0].edges[edgeIndex0].edgeNumber;

            IndexType cellIndex1 = contactEdgesLocation[1].edges[edgeIndex1].cellIndex;
            IndexType edgeNumber1 = contactEdgesLocation[1].edges[edgeIndex1].edgeNumber;

            additionalCellInfos[cellIndex0].neighbouringEdges[edgeNumber0].correspondingCellIndex = cellIndex1;
            additionalCellInfos[cellIndex0].neighbouringEdges[edgeNumber0].correspondingEdgeNumber = edgeNumber1;
            additionalCellInfos[cellIndex0].neighbouringEdges[edgeNumber0].interactionType = IndexType(contactTypeIndex);

            additionalCellInfos[cellIndex1].neighbouringEdges[edgeNumber1].correspondingCellIndex = cellIndex0;
            additionalCellInfos[cellIndex1].neighbouringEdges[edgeNumber1].correspondingEdgeNumber = edgeNumber0;
            additionalCellInfos[cellIndex1].neighbouringEdges[edgeNumber1].interactionType = IndexType(contactTypeIndex);
          }
        }
      }
    }
    offset += contactEdgesCount[contactTypeIndex];
  }
}

void GeomMesh<Space2>::BuildBoundariesInfo(BoundaryEdge* boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount)
{
  IndexType offset = 0;
  for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; ++boundaryTypeIndex)
  {    
    for (IndexType boundaryEdgeIndex = 0; boundaryEdgeIndex < boundaryEdgesCount[boundaryTypeIndex]; ++boundaryEdgeIndex)
    {
      BoundaryEdge boundaryEdge = boundaryEdges[offset + boundaryEdgeIndex];  
      IndexType nodeIndex1 = boundaryEdge.nodeIndices[0];
      IndexType nodeIndex2 = boundaryEdge.nodeIndices[1];
      if (nodeIndex1 != IndexType(-1) && nodeIndex2 != IndexType(-1)) 
      {
        EdgeLocationPair edgeLocationPair = GetEdgeLocation(nodeIndex1, nodeIndex2);
        IndexType cellIndex1 = edgeLocationPair.edges[0].cellIndex;
        IndexType cellIndex2 = edgeLocationPair.edges[1].cellIndex;
        IndexType edgeNumber1 = edgeLocationPair.edges[0].edgeNumber;
        IndexType edgeNumber2 = edgeLocationPair.edges[1].edgeNumber;

        if (!edgeLocationPair.edges[0].IsNull())
        {                      
          additionalCellInfos[cellIndex1].neighbouringEdges[edgeNumber1].correspondingCellIndex = IndexType(-1);
          additionalCellInfos[cellIndex1].neighbouringEdges[edgeNumber1].correspondingEdgeNumber = IndexType(-1);
          additionalCellInfos[cellIndex1].neighbouringEdges[edgeNumber1].interactionType = boundaryTypeIndex;
        }

        if (!edgeLocationPair.edges[1].IsNull())
        {
          additionalCellInfos[cellIndex2].neighbouringEdges[edgeNumber2].correspondingCellIndex = IndexType(-1);
          additionalCellInfos[cellIndex2].neighbouringEdges[edgeNumber2].correspondingEdgeNumber = IndexType(-1);
          additionalCellInfos[cellIndex2].neighbouringEdges[edgeNumber2].interactionType = boundaryTypeIndex;
        }
      }
    }
    offset += boundaryEdgesCount[boundaryTypeIndex];
  }
}

bool GeomMesh<Space2>::FindCommonEdge(IndexType srcCellIndex, IndexType dstCellIndex, 
  IndexType &_srcEdgeNumber, IndexType &_dstEdgeNumber)
{
  IndexType srcCellNodes[3];
  IndexType dstCellNodes[3];

  GetFixedCellIndices(srcCellIndex, srcCellNodes);
  GetFixedCellIndices(dstCellIndex, dstCellNodes);

  for(IndexType srcEdgeNumber = 0; srcEdgeNumber < 3; srcEdgeNumber++)
  {
    IndexType srcEdgeNodes[3];
    GetCellEdgeNodes(srcCellNodes, srcEdgeNumber, srcEdgeNodes);

    for(IndexType dstEdgeNumber = 0; dstEdgeNumber < 3; dstEdgeNumber++)
    {
      IndexType dstEdgeNodes[3];
      GetCellEdgeNodes(dstCellNodes, dstEdgeNumber, dstEdgeNodes);

      if(srcEdgeNodes[0] == dstEdgeNodes[1] && srcEdgeNodes[1] == dstEdgeNodes[0])
      {
        _srcEdgeNumber = srcEdgeNumber;
        _dstEdgeNumber = dstEdgeNumber;
        return true;
      }
    }
  }
  return false;
}

void GeomMesh<Space2>::BuildCellAdditionalTopology(IndexType *internalContactTypes)
{
  additionalCellInfos.resize(cells.size());

  for (IndexType edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++)
  {
    IndexType srcEdgeNumber;
    IndexType dstEdgeNumber;

    if (edgeTopologyInfos[edgeIndex].incidentCellsCount == 2)
    {
      IndexType contactingCells[2];
      for (IndexType contactingCellIndex = 0; contactingCellIndex < 2; contactingCellIndex++)
      {
        contactingCells[contactingCellIndex] = edgeTopologyInfos[edgeIndex].incidentCells[contactingCellIndex];
      }

      for (IndexType cellNumber = 0; cellNumber < 2; cellNumber++)
      {
        IndexType srcCellIndex = contactingCells[cellNumber];
        IndexType dstCellIndex = contactingCells[(cellNumber + 1) % 2];

        if (FindCommonEdge(srcCellIndex, dstCellIndex, srcEdgeNumber, dstEdgeNumber))
        {
          additionalCellInfos[srcCellIndex].neighbouringEdges[srcEdgeNumber].correspondingCellIndex = dstCellIndex;
          additionalCellInfos[srcCellIndex].neighbouringEdges[srcEdgeNumber].correspondingEdgeNumber = dstEdgeNumber;
          if(internalContactTypes)
            additionalCellInfos[srcCellIndex].neighbouringEdges[srcEdgeNumber].interactionType = 
              std::max(internalContactTypes[srcCellIndex], internalContactTypes[dstCellIndex]); //priority is given to a larger type
          else
            additionalCellInfos[srcCellIndex].neighbouringEdges[srcEdgeNumber].interactionType = IndexType(0);
        } else
        {
          assert(0); //topology error
        }
      }
    }
  }
}

void GeomMesh<Space2>::GetCellFaceNodes(IndexType cellIndex, IndexType faceNumber, IndexType* faceNodes) const
{
  GetCellEdgeNodes(cellIndex, faceNumber, faceNodes);
}

void GeomMesh<Space2>::GetCellEdgeVertices(const IndexType cellIndex, IndexType edgeNumber, Vector* edgeVertices) const
{
  Vector cellVertices[Space2::NodesPerCell];
  GetCellVertices(cellIndex, cellVertices);
  edgeVertices[0] = cellVertices[edgeNumber];
  edgeVertices[1] = cellVertices[(edgeNumber + 1) % 3];
}

Space2::Vector GeomMesh<Space2>::GetCellEdgeMiddle(const IndexType cellIndex, IndexType edgeNumber) const
{
  Vector edgeVertices[2];
  GetCellEdgeVertices(cellIndex, edgeNumber, edgeVertices);
  return (edgeVertices[0] + edgeVertices[1]) * Scalar(0.5);
}

void GeomMesh<Space2>::GetGhostCellVertices(IndexType cellIndex, IndexType boundaryEdgeNumber, Vector* ghostCellVertices) const
{
  Vector points[Space::NodesPerCell];
  GetCellVertices(cellIndex, points);

  for (IndexType vertexIndex = 0; vertexIndex < Space::NodesPerCell; ++vertexIndex)
  {
    ghostCellVertices[vertexIndex] = points[vertexIndex];
  }

  Vector fixedEdge = points[boundaryEdgeNumber] - points[(boundaryEdgeNumber + 1) % 3];
  Vector edgeNormal = fixedEdge.GetNorm().GetPerpendicular();

  IndexType mirroringNodeNumber = 3 - boundaryEdgeNumber - (boundaryEdgeNumber + 1) % 3;
  
  //Vector2 reflectedPoint = globalPoint - edgeNormal * (edgeNormal * (globalPoint - edgePoint)) * Scalar(2.0);

  ghostCellVertices[mirroringNodeNumber] = points[mirroringNodeNumber] - 
                                           edgeNormal * (edgeNormal * (points[mirroringNodeNumber] - points[boundaryEdgeNumber])) * Scalar(2.0);

  //std::swap(ghostCellVertices[boundaryEdgeNumber], ghostCellVertices[(boundaryEdgeNumber + 1) % 3]);
}

Space2::Vector GeomMesh<Space2>::GetEdgeExternalNormal(IndexType cellIndex, IndexType edgeNumber) const
{
  Vector points[Space2::NodesPerCell];
  GetCellVertices(cellIndex, points);
  Vector edge = points[(edgeNumber + 1) % 3] - points[edgeNumber];
  return -edge.GetPerpendicular().GetNorm();
}

Space2::Vector GeomMesh<Space2>::GetFaceExternalNormal(IndexType cellIndex, IndexType faceNumber) const
{
  return GetEdgeExternalNormal(cellIndex, faceNumber);
}

GeomMesh<Space2>::EdgeLocationPair GeomMesh<Space2>::BuildEdgeLocation(const EdgePairIndices& edgePairIndices)
{
  EdgeLocationPair edgeLocationPair;
  IndexType edgesFoundCount = 0;

  for (IndexType edgePairIndex = 0; edgePairIndex < 2; ++edgePairIndex)
  {
    IndexType nodeIndex0 = edgePairIndices.edges[edgePairIndex].nodeIndices[0];
    IndexType nodeIndex1 = edgePairIndices.edges[edgePairIndex].nodeIndices[1];
    EdgeLocationPair edgesLocation = BuildEdgeLocation(nodeIndex0, nodeIndex1);

    for (IndexType edgeIndex = 0; edgeIndex < 2; ++edgeIndex)
    {
      if (!edgesLocation.edges[edgeIndex].IsNull())
      {
        bool found = false;
        for (IndexType edgesFoundIndex = 0; edgesFoundIndex < edgesFoundCount; ++edgesFoundIndex)
        {
          if (edgeLocationPair.edges[edgesFoundIndex] == edgesLocation.edges[edgeIndex]) found = true;
        }
        if (!found) edgeLocationPair.edges[edgesFoundCount++] = edgesLocation.edges[edgeIndex];
      }
    }
  }
  assert(edgesFoundCount == 2);
  return edgeLocationPair;
}

// to replace
bool GeomMesh<Space2>::IsBoundaryCell(IndexType cellIndex) const
{
  for (IndexType edgeNumber = 0; edgeNumber < Space2::FacesPerCell; ++edgeNumber)
  {
    if (additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex == IndexType(-1))
    {
      return true;
    }
  }
  return false;
}

bool GeomMesh<Space2>::IsBoundaryNode(IndexType nodeIndex) const
{
  for (IndexType incidentCellNumber = 0; incidentCellNumber < nodeTopologyInfos[nodeIndex].incidentCellsCount; ++incidentCellNumber)
  {
    IndexType incidentCellIndex = nodeTopologyInfos[nodeIndex].incidentCells[incidentCellNumber];
    IndexType nodeIndices[3];
    GetFixedCellIndices(incidentCellIndex, nodeIndices);

    for (IndexType incidentNodeNumber = 0; incidentNodeNumber < 3; ++incidentNodeNumber)
    {
      if (nodeIndices[incidentNodeNumber] != nodeIndex)
      {
        EdgeLocationPair locations = GetEdgeLocation(nodeIndex, nodeIndices[incidentNodeNumber]);
        for (IndexType pairIndex = 0; pairIndex < 2; ++pairIndex)
        {
          if (!locations.edges[pairIndex].IsNull())
          {
            EdgeLocation edge = locations.edges[pairIndex];
            if (additionalCellInfos[edge.cellIndex].neighbouringEdges[edge.edgeNumber].correspondingCellIndex == IndexType(-1))
            {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}
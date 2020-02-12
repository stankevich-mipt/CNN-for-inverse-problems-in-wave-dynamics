#include <stack>
#include <set>
#include <algorithm>

// Node
template <typename Space>
GeomMeshCommon<Space>::Node::Node()
{
}

template <typename Space>
GeomMeshCommon<Space>::Node::Node(const Vector& pos): pos(pos)
{
}

// NodeTopology
template <typename Space>
GeomMeshCommon<Space>::NodeTopologyInfo::NodeTopologyInfo():
  incidentCellsCount(0),
  incidentCells(0)
{
}

// Cell
template <typename Space>
GeomMeshCommon<Space>::Cell::Cell()
{
  std::fill_n(incidentNodes, Space::NodesPerCell, IndexType(-1));
}

template <typename Space>
GeomMeshCommon<Space>::Cell::Cell(IndexType nodeIndices[Space::NodesPerCell])
{
  std::copy(nodeIndices, nodeIndices + Space::NodesPerCell, incidentNodes);
  std::sort(incidentNodes, incidentNodes + Space::NodesPerCell);
}

template <typename Space>
bool GeomMeshCommon<Space>::Cell::operator==(const Cell& other) const
{
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    if (incidentNodes[nodeNumber] != other.incidentNodes[nodeNumber]) return false;
  }
  return true;
}

template <typename Space>
bool GeomMeshCommon<Space>::Cell::operator<(const Cell& other) const
{
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    if (incidentNodes[nodeNumber] == other.incidentNodes[nodeNumber]) continue;
    return incidentNodes[nodeNumber] < other.incidentNodes[nodeNumber];
  }
  return false;
}

template <typename Space>
bool GeomMeshCommon<Space>::Cell::operator!=(const Cell& other) const
{
  return !this->operator==(other);
}

template <typename Space>
void GeomMeshCommon<Space>::Cell::Copy(const Cell& other)
{
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    incidentNodes[nodeNumber] = other.incidentNodes[nodeNumber];
  }
}

template <typename Space>
GeomMeshCommon<Space>::Cell::Cell(const Cell& other)
{
  Copy(other);
}

template <typename Space>
typename GeomMeshCommon<Space>::Cell& GeomMeshCommon<Space>::Cell::operator=(const Cell& other)
{
  Copy(other);
  return *this;
}

// Edge
template <typename Space>
GeomMeshCommon<Space>::Edge::Edge()
{
  std::fill_n(incidentNodes, Space::NodesPerEdge, IndexType(-1));
}

template <typename Space>
GeomMeshCommon<Space>::Edge::Edge(IndexType nodeIndices[Space::NodesPerEdge])
{
  std::copy(nodeIndices, nodeIndices + Space::NodesPerEdge, incidentNodes);
  std::sort(incidentNodes, incidentNodes + Space::NodesPerEdge);
}

template <typename Space>
GeomMeshCommon<Space>::Edge::Edge(IndexType nodeIndex0, IndexType nodeIndex1)
{
  incidentNodes[0] = nodeIndex0;
  incidentNodes[1] = nodeIndex1;
  std::sort(incidentNodes, incidentNodes + Space::NodesPerEdge);
}

template <typename Space>
bool GeomMeshCommon<Space>::Edge::operator==(const Edge& other) const
{
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
  {
    if (incidentNodes[nodeNumber] != other.incidentNodes[nodeNumber]) return false;
  }
  return true;
}

template <typename Space>
bool GeomMeshCommon<Space>::Edge::operator<(const Edge& other) const
{
  return (incidentNodes[0]  < other.incidentNodes[0]) || 
         (incidentNodes[0] == other.incidentNodes[0] && incidentNodes[1] < other.incidentNodes[1]);
}

template <typename Space>
bool GeomMeshCommon<Space>::Edge::operator!=(const Edge& other) const
{
  return !this->operator==(other);
}

template <typename Space>
void GeomMeshCommon<Space>::Edge::Copy(const Edge& other)
{
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
  {
    incidentNodes[nodeNumber] = other.incidentNodes[nodeNumber];
  }
}

template <typename Space>
GeomMeshCommon<Space>::Edge::Edge(const Edge& other)
{
  Copy(other);
}

template <typename Space>
typename GeomMeshCommon<Space>::Edge& GeomMeshCommon<Space>::Edge::operator=(const Edge& other)
{
  Copy(other);
  return *this;
}

// Edge topology
template <typename Space>
GeomMeshCommon<Space>::EdgeTopologyInfo::EdgeTopologyInfo():
  incidentCellsCount(0),
  incidentCells(0)
{
}

template <typename Space>
GeomMeshCommon<Space>::EdgeIndices::EdgeIndices()
{
  std::fill_n(nodeIndices, Space::NodesPerEdge, IndexType(-1));
}

template <typename Space>
GeomMeshCommon<Space>::EdgeIndices::EdgeIndices(IndexType nodeIndex0, IndexType nodeIndex1)
{
  nodeIndices[0] = nodeIndex0;
  nodeIndices[1] = nodeIndex1;
  std::sort(nodeIndices, nodeIndices + Space::NodesPerEdge);
}

template <typename Space>
bool GeomMeshCommon<Space>::EdgeIndices::operator<(const EdgeIndices& other) const
{
  return nodeIndices[0] < other.nodeIndices[0] || 
          (nodeIndices[0] == other.nodeIndices[0] && nodeIndices[1] < other.nodeIndices[1]);
}

template <typename Space>
bool GeomMeshCommon<Space>::EdgeIndices::operator==(const EdgeIndices& other) const
{
  return nodeIndices[0] == other.nodeIndices[0] && nodeIndices[1] == other.nodeIndices[1];
}

template <typename Space>
bool GeomMeshCommon<Space>::EdgeIndices::operator!=(const EdgeIndices& other) const
{
  return !this->operator==(other);
}

template <typename Space>
void GeomMeshCommon<Space>::EdgeIndices::Copy(const EdgeIndices& other)
{
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
  {
    nodeIndices[nodeNumber] = other.nodeIndices[nodeNumber];
  }
}

template <typename Space>
GeomMeshCommon<Space>::EdgeIndices::EdgeIndices(const EdgeIndices& other)
{
  Copy(other);
}

template <typename Space>
typename GeomMeshCommon<Space>::EdgeIndices& GeomMeshCommon<Space>::EdgeIndices::operator=(const EdgeIndices& other)
{
  Copy(other);
  return *this;
}

template <typename Space>
bool GeomMeshCommon<Space>::EdgePairIndices::IsEdgesCoincident() const
{
  return edges[0] == edges[1];
}

template <typename Space>
bool GeomMeshCommon<Space>::EdgePairIndices::operator==(const EdgePairIndices& other) const
{
  return edges[0] == other.edges[0] && edges[1] == other.edges[1];
}

template <>
void GeomMeshCommon<Space2>::GetCellEdgeNodes(const IndexType* cellIncidentNodes, IndexType edgeNumber, IndexType* edgeNodeIndices) const
{
  edgeNodeIndices[0] = cellIncidentNodes[edgeNumber];
  edgeNodeIndices[1] = cellIncidentNodes[(edgeNumber + 1) % 3];
}

template <>
void GeomMeshCommon<Space2>::GetCellIndices(IndexType cellIndex, IndexType* incidentNodes) const
{
  std::copy(cells[cellIndex].incidentNodes, cells[cellIndex].incidentNodes + Space::NodesPerCell, incidentNodes);
}

template <>
void GeomMeshCommon<Space2>::GetFixedCellIndices(IndexType cellIndex, IndexType* incidentNodes) const
{
  GetCellIndices(cellIndex, incidentNodes);
  if(((nodes[incidentNodes[1]].pos - nodes[incidentNodes[0]].pos) ^ (nodes[incidentNodes[2]].pos - nodes[incidentNodes[0]].pos)) < 0)
  {
    std::swap(incidentNodes[0], incidentNodes[1]);
  }
}

template <>
void GeomMeshCommon<Space3>::GetCellIndices(IndexType cellIndex, IndexType* incidentNodes) const
{
  std::copy(cells[cellIndex].incidentNodes, cells[cellIndex].incidentNodes + Space::NodesPerCell, incidentNodes);
}

template <>
void GeomMeshCommon<Space3>::GetFixedCellIndices(IndexType cellIndex, IndexType* incidentNodes) const
{
  GetCellIndices(cellIndex, incidentNodes);
  if(MixedProduct(nodes[incidentNodes[1]].pos - nodes[incidentNodes[0]].pos,
                  nodes[incidentNodes[2]].pos - nodes[incidentNodes[0]].pos,
                  nodes[incidentNodes[3]].pos - nodes[incidentNodes[0]].pos) < 0)
  {
    std::swap(incidentNodes[0], incidentNodes[1]);
  }
}

template <>
void GeomMeshCommon<Space2>::GetCellEdgeNodes(const IndexType cellIndex, IndexType edgeNumber, IndexType* edgeNodes) const
{
  IndexType cellIncidentNodes[Space::NodesPerCell];
  GetFixedCellIndices(cellIndex, cellIncidentNodes);
  return GetCellEdgeNodes(cellIncidentNodes, edgeNumber, edgeNodes);
}

template <>
void GeomMeshCommon<Space3>::GetCellEdgeNodes(const IndexType* cellIncidentNodes, IndexType edgeNumber, IndexType* edgeNodes) const
{
  assert(0);
  // TODO
}

template <>
void GeomMeshCommon<Space3>::GetCellEdgeNodes(const IndexType cellIndex, IndexType edgeNumber, IndexType* edgeNodes) const
{
  assert(0);
  // TODO
}

template <>
void GeomMeshCommon<Space2>::LoadGeom(Vector* vertices, IndexType* indices, IndexType verticesCount, IndexType cellsCount)
{
  // cells
  cells.reserve(cellsCount);
  for(IndexType cellIndex = 0; cellIndex < cellsCount; cellIndex++)
  {
    IndexType cellNodeIndices[Space::NodesPerCell];
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      cellNodeIndices[nodeNumber] = indices[cellIndex * Space::NodesPerCell + nodeNumber];
    }
    cells.push_back(Cell(cellNodeIndices));
  }
  std::sort(cells.begin(), cells.end());

  // nodes
  nodes.reserve(verticesCount);
  for(IndexType nodeIndex = 0; nodeIndex < verticesCount; nodeIndex++)
  {
    Node newbie;
    newbie.pos = vertices[nodeIndex];
    nodes.push_back(newbie);
  }

  // edges
  for(IndexType cellIndex = 0; cellIndex < cellsCount; cellIndex++)
  {
    for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
    {
      IndexType edgeNodeIndices[Space::NodesPerEdge];
      GetCellEdgeNodes(cells[cellIndex].incidentNodes, edgeNumber, edgeNodeIndices);
      edges.push_back(Edge(edgeNodeIndices));
    }
  }
  std::sort(edges.begin(), edges.end());
  RemoveDuplicates(&edges);
}

template <>
void GeomMeshCommon<Space3>::BuildTopologyInfos()
{
}

template <>
void GeomMeshCommon<Space2>::BuildNodeTopologyInfos()
{
  nodeTopologyInfos.resize(nodes.size());

  IndexType totalNodeIncidentCells = 0;

  for(IndexType phase = 0; phase < 2; phase++)
  {
    if (phase == 1)
    {
      nodeIncidentCellsPool.resize(totalNodeIncidentCells);
      IndexType offset = 0;
      for (IndexType nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex)
      {
        nodeTopologyInfos[nodeIndex].incidentCells = nodeTopologyInfos[nodeIndex].incidentCellsCount ? &(nodeIncidentCellsPool[offset]) : 0;
        offset += nodeTopologyInfos[nodeIndex].incidentCellsCount;
        nodeTopologyInfos[nodeIndex].incidentCellsCount = 0;
      }
      assert(offset == totalNodeIncidentCells);
    }
 
    for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex) 
    {
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber) 
      {
        IndexType nodeIndex = cells[cellIndex].incidentNodes[nodeNumber];
        if(phase == 0)
        {
          nodeTopologyInfos[nodeIndex].incidentCellsCount++;
          totalNodeIncidentCells++;
        } else
        {
          nodeTopologyInfos[nodeIndex].incidentCells[nodeTopologyInfos[nodeIndex].incidentCellsCount++] = cellIndex;
        }
      }
    }
  }
}

template <>
void GeomMeshCommon<Space2>::BuildEdgeTopologyInfos()
{
  edgeTopologyInfos.resize(edges.size());
  IndexType totalEdgeIncidentCells = 0;

  for (IndexType phase = 0; phase < 2; ++phase)
  {
    if (phase == 1)
    {
      edgeIncidentCellsPool.resize(totalEdgeIncidentCells);
      IndexType offset = 0;
      for (IndexType edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex)
      {
        edgeTopologyInfos[edgeIndex].incidentCells = edgeTopologyInfos[edgeIndex].incidentCellsCount ? &(edgeIncidentCellsPool[offset]) : 0;
        offset += edgeTopologyInfos[edgeIndex].incidentCellsCount;
        edgeTopologyInfos[edgeIndex].incidentCellsCount = 0;
      }
      assert(offset == totalEdgeIncidentCells);
    }

    for (IndexType edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex)
    {
      IndexType baseNode = edges[edgeIndex].incidentNodes[0];
      for (IndexType cellNumber = 0; cellNumber < nodeTopologyInfos[baseNode].incidentCellsCount; ++cellNumber)
      {
        IndexType incidentCellIndex = nodeTopologyInfos[baseNode].incidentCells[cellNumber];
        for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; ++edgeNumber)
        {
          IndexType edgeNodeIndices[Space::NodesPerEdge];
          GetCellEdgeNodes(cells[incidentCellIndex].incidentNodes, edgeNumber, edgeNodeIndices);
          if (edges[edgeIndex] == Edge(edgeNodeIndices))
          {
            if (phase == 0)
            {
              edgeTopologyInfos[edgeIndex].incidentCellsCount++;
              totalEdgeIncidentCells++;
            } else
            {
              edgeTopologyInfos[edgeIndex].incidentCells[edgeTopologyInfos[edgeIndex].incidentCellsCount++] = incidentCellIndex;
            }
          }
        }
      }
    }
  }
}

template <>
void GeomMeshCommon<Space2>::BuildCellTopologyInfos()
{
  cellTopologyInfos.resize(cells.size());
  std::vector<IndexType> cellIncidentEdgesCount(cells.size());

  for (IndexType edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex)
  {
    for (IndexType cellNumber = 0; cellNumber < edgeTopologyInfos[edgeIndex].incidentCellsCount; ++cellNumber) 
    {
      IndexType incidentCellIndex = edgeTopologyInfos[edgeIndex].incidentCells[cellNumber];
      cellTopologyInfos[incidentCellIndex].incidentEdges[cellIncidentEdgesCount[incidentCellIndex]++] = edgeIndex;
    }
  }
}

template <>
void GeomMeshCommon<Space2>::BuildTopologyInfos()
{
  BuildNodeTopologyInfos();
  BuildEdgeTopologyInfos();
  BuildCellTopologyInfos();
}

template <>
void GeomMeshCommon<Space3>::LoadGeom(Vector* vertexPositions, IndexType* cellIndices, IndexType verticesCount, IndexType indicesCount)
{
  TopologyReconstructor<IndexType> rec(cellIndices, indicesCount);

  IndexType nodesCount = rec.GetNodesCount();
  IndexType edgesCount = rec.GetEdgesCount();
  IndexType facesCount = rec.GetFacesCount();
  IndexType cellsCount = rec.GetCellsCount();

  nodes.resize(nodesCount);
  edges.resize(edgesCount);
  faces.resize(facesCount);
  cells.resize(cellsCount);

  //nodes
  IndexType totalNodeIncidentCells = rec.GetTotalNodeIncidentCells();
  nodeIncidentCellsPool.resize(totalNodeIncidentCells);

  IndexType nodePoolOffset = 0;
  nodeTopologyInfos.resize(nodesCount);

  for(IndexType i = 0; i < nodesCount; i++)
  {
    nodes[i].pos = vertexPositions[i];

    nodeTopologyInfos[i].incidentCells = &(nodeIncidentCellsPool[nodePoolOffset]);
    rec.GetNodeIncidentCells(i, nodeTopologyInfos[i].incidentCells);
    nodeTopologyInfos[i].incidentCellsCount = rec.GetNodeIncidentCellsCount(i);
    nodePoolOffset += rec.GetNodeIncidentCellsCount(i);
  }

  //edges
  IndexType totalEdgeIncidentCells = rec.GetTotalEdgeIncidentCells();
  edgeIncidentCellsPool.resize(totalEdgeIncidentCells);
  edgeTopologyInfos.resize(edgesCount);

  IndexType edgePoolOffset = 0;

  for(IndexType i = 0; i < edgesCount; i++)
  {
    rec.GetEdgeIncidentNodes(i, edges[i].incidentNodes);

    edgeTopologyInfos[i].incidentCells = &(edgeIncidentCellsPool[edgePoolOffset]);
    rec.GetEdgeIncidentCells(i, edgeTopologyInfos[i].incidentCells);
    edgeTopologyInfos[i].incidentCellsCount = rec.GetEdgeIncidentCellsCount(i);
    edgePoolOffset += rec.GetEdgeIncidentCellsCount(i);
  }

  //faces
  IndexType totalFaceIncidentCells = rec.GetTotalFaceIncidentCells();
  faceIncidentCellsPool.resize(totalFaceIncidentCells);
  faceTopologyInfos.resize(facesCount);

  IndexType facePoolOffset = 0;

  for(IndexType i = 0; i < facesCount; i++)
  {
    rec.GetFaceIncidentNodes(i, faces[i].incidentNodes);

    faceTopologyInfos[i].incidentCells = &(faceIncidentCellsPool[facePoolOffset]);
    rec.GetFaceIncidentCells(i, faceTopologyInfos[i].incidentCells);
    faceTopologyInfos[i].incidentCellsCount = rec.GetFaceIncidentCellsCount(i);
    facePoolOffset += rec.GetFaceIncidentCellsCount(i);
  }

  //cells
  cellTopologyInfos.resize(cellsCount);
  for(IndexType i = 0; i < cellsCount; i++)
  {
    rec.GetCellIncidentNodes(i, cells[i].incidentNodes);
    rec.GetCellIncidentEdges(i, cellTopologyInfos[i].incidentEdges);
    rec.GetCellIncidentFaces(i, cellTopologyInfos[i].incidentFaces);
  }
}

template <typename Space>
void GeomMeshCommon<Space>::GetCellVertices(IndexType cellIndex, Vector* points) const
{
  IndexType cellIndices[Space::NodesPerCell];
  GetFixedCellIndices(cellIndex, cellIndices);

  GetCellVertices(cellIndices, points);
}

template <typename Space>
void GeomMeshCommon<Space>::GetCellVertices(const IndexType* cellIncidentNodes, Vector* points) const
{
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; nodeNumber++)
  {
    points[nodeNumber] = nodes[cellIncidentNodes[nodeNumber]].pos;
  }
}

template <typename Space>
typename Space::Vector GeomMeshCommon<Space>::GetMassCenter(IndexType cellIndex) const
{
  Vector center = Vector::zero();
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    IndexType nodeIndex = cells[cellIndex].incidentNodes[nodeNumber];
    center += nodes[nodeIndex].pos;
  }
  return center / Scalar(Space::NodesPerCell);
}

// mass center of all nodes
template <typename Space>
typename Space::Vector GeomMeshCommon<Space>::GetMassCenter() const
{
  Vector center = Vector::zero();

  for (IndexType nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex)
  {
    center += nodes[nodeIndex].pos;
  }
  return center / Scalar(nodes.size());
}

template <>
Space2::Scalar GeomMeshCommon<Space2>::GetMinHeight(IndexType cellIndex) const
{
  Vector points[Space2::NodesPerCell];
  GetCellVertices(cellIndex, points);

  Scalar minHeight = std::numeric_limits<Scalar>::max();

  for (IndexType topVertexNumber = 0; topVertexNumber < Space2::NodesPerCell; ++topVertexNumber)
  {
    Vector edges[Space2::Dimension];
    IndexType edgesCount = 0;
    for (IndexType vertexNumber = 0; vertexNumber < Space2::NodesPerCell; ++vertexNumber)
    {
      if (topVertexNumber != vertexNumber)
      {
        edges[edgesCount++] = points[vertexNumber] - points[topVertexNumber];
      }
    }
    minHeight = std::min(minHeight, TriHeight(edges));
  }
  return minHeight;
}

template <>
Space3::Scalar GeomMeshCommon<Space3>::GetMinHeight(IndexType cellIndex) const
{
  Vector points[Space3::NodesPerCell];
  GetCellVertices(cellIndex, points);

  Scalar minHeight = std::numeric_limits<Scalar>::max();

  for (IndexType topVertexNumber = 0; topVertexNumber < Space3::NodesPerCell; ++topVertexNumber)
  {
    Vector edges[Space3::Dimension];
    IndexType edgesCount = 0;
    for (IndexType vertexNumber = 0; vertexNumber < Space3::NodesPerCell; ++vertexNumber)
    {
      if (topVertexNumber != vertexNumber)
      {
        edges[edgesCount++] = points[vertexNumber] - points[topVertexNumber];
      }
    }
    minHeight = std::min(minHeight, TetraHeight(edges));
  }
  return minHeight;
}

template <typename Space>
typename Space::Scalar GeomMeshCommon<Space>::GetMinHeight() const
{
  Scalar minHeight = std::numeric_limits<Scalar>::max();
  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    minHeight = Min(minHeight, GetMinHeight(cellIndex));
  }
  return minHeight;
}

template <typename Space>
typename Space::AABB GeomMeshCommon<Space>::GetCellAABB(IndexType cellIndex) const
{
  Vector vertices[Space::NodesPerCell];
  GetCellVertices(cellIndex, vertices);

  AABB aabb;
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    aabb.Expand(vertices[nodeNumber]);
  }
  return aabb;
}

template <typename Space>
void GeomMeshCommon<Space>::UpdateAABBTree()
{
  for (IndexType cellIndex = 0; cellIndex < cells.size(); cellIndex++)
  if (IsCellInAABBTree(cellIndex))
  {
    aabbTree.UpdateNode(treeNodeCellIndices[cellIndex], GetCellAABB(cellIndex));
  }
}

template <typename Space>
void GeomMeshCommon<Space>::BuildAABBTree()
{
  printf("Building AABB tree for triangles\n");
  treeNodeCellIndices.reserve(cells.size());
  for (IndexType cellIndex = 0; cellIndex < cells.size(); cellIndex++)
  {
    treeNodeCellIndices.push_back(aabbTree.InsertNode(GetCellAABB(cellIndex)));
    aabbTree.SetUserData(treeNodeCellIndices[cellIndex], cellIndex);
  }
}

template <typename Space>
typename Space::AABB GeomMeshCommon<Space>::GetAABB() const
{
  AABB aabb;
  for (IndexType nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex)
  {
    aabb.Expand(nodes[nodeIndex].pos);
  }
  return aabb;
}

template <typename Space>
void GeomMeshCommon<Space>::RemoveCellFromAABBTree(IndexType cellIndex)
{ 
  if (!treeNodeCellIndices.empty())
  {
    IndexType treeNodeCellIndex = treeNodeCellIndices[cellIndex];
    if (treeNodeCellIndex != IndexType(-1))
    {
      treeNodeCellIndices[cellIndex] = IndexType(-1);
      aabbTree.RemoveNode(treeNodeCellIndex);
    }
  }
}

template <typename Space>
void GeomMeshCommon<Space>::AddToAABBTree(IndexType cellIndex)
{
  if (!treeNodeCellIndices.empty())
  {
    treeNodeCellIndices[cellIndex] = aabbTree.InsertNode(GetCellAABB(cellIndex));
    aabbTree.SetUserData(treeNodeCellIndices[cellIndex], cellIndex);
  }
}

template <typename Space>
bool GeomMeshCommon<Space>::IsCellInAABBTree(IndexType cellIndex) const
{
  return !treeNodeCellIndices.empty() && treeNodeCellIndices[cellIndex] != IndexType(-1);
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::GetIncidentCellsCount(IndexType nodeIndex) const
{
  return nodeTopologyInfos[nodeIndex].incidentCellsCount;
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::GetIncidentCellIndex(IndexType nodeIndex, IndexType cellNumber) const
{
  return nodeTopologyInfos[nodeIndex].incidentCells[cellNumber];
}

template <>
Space2::IndexType GeomMeshCommon<Space2>::GetCorrespondingCellIndex(IndexType cellIndex, IndexType edgeNumber) const
{
  return additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex;
}

template <>
Space3::IndexType GeomMeshCommon<Space3>::GetCorrespondingCellIndex(IndexType cellIndex, IndexType faceNumber) const
{
  return additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
}

template <>
Space2::IndexType GeomMeshCommon<Space2>::GetCorrespondingFaceNumber(IndexType cellIndex, IndexType edgeNumber) const
{
  return additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber;
}

template <>
Space3::IndexType GeomMeshCommon<Space3>::GetCorrespondingFaceNumber(IndexType cellIndex, IndexType faceNumber) const
{
  return additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber;
}


template <>
Space2::IndexType GeomMeshCommon<Space2>::GetInteractionType(IndexType cellIndex, IndexType edgeNumber) const
{
  return additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;
}

template <>
Space3::IndexType GeomMeshCommon<Space3>::GetInteractionType(IndexType cellIndex, IndexType faceNumber) const
{
  return additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;
}

template <>
void GeomMeshCommon<Space3>::BuildNodeTopologyInfos()
{
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::GetCellIndex(IndexType nodeIndices[Space::NodesPerCell]) const
{
  IndexType orderedIndices[Space::NodesPerCell];
  std::copy(nodeIndices, nodeIndices + Space::NodesPerCell, orderedIndices);
  std::sort(orderedIndices, orderedIndices + Space::NodesPerCell);

  IndexType baseNode = orderedIndices[0]; //may be any of given nodes
  for(IndexType incidentCellIndex = 0; incidentCellIndex < nodeTopologyInfos[baseNode].incidentCellsCount; incidentCellIndex++)
  {
    IndexType cellIndex = nodeTopologyInfos[baseNode].incidentCells[incidentCellIndex];
    bool found = true;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (cells[cellIndex].incidentNodes[nodeNumber] != orderedIndices[nodeNumber])
      {
        found = false;
        break;
      }
    }
    if (found)
    {
      return cellIndex;
    }
  }
  return IndexType(-1);
}

template <typename Space>
typename GeomMeshCommon<Space>::Node* GeomMeshCommon<Space>::GetNode(IndexType index)
{
  return &(nodes[index]);
}

template <typename Space>
const typename GeomMeshCommon<Space>::Node* GeomMeshCommon<Space>::GetNode(IndexType index) const
{
  return &(nodes[index]);
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::GetNodesCount()
{
  return IndexType(nodes.size());
}

template <typename Space>
typename GeomMeshCommon<Space>::Edge* GeomMeshCommon<Space>::GetEdge(IndexType index)
{
  return &(edges[index]);
}

template <typename Space>
const typename GeomMeshCommon<Space>::Edge* GeomMeshCommon<Space>::GetEdge(IndexType index) const
{
  return &(edges[index]);
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::GetEdgesCount()
{
  return IndexType(edges.size());
}

template <typename Space>
typename GeomMeshCommon<Space>::Cell* GeomMeshCommon<Space>::GetCell(IndexType index)
{
  return &(cells[index]);
}

template <typename Space>
const typename GeomMeshCommon<Space>::Cell* GeomMeshCommon<Space>::GetCell(IndexType index) const
{
  return &(cells[index]);
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::GetCellsCount()
{
  return IndexType(cells.size());
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::
  GetEdgeIndex(IndexType nodeIndices[Space::NodesPerEdge]) const
{
  std::sort(nodeIndices, nodeIndices + Space::NodesPerEdge);
  IndexType baseNode = nodeIndices[0]; //may be any of given nodes

  for(IndexType incidentCellIndex = 0; incidentCellIndex < nodeTopologyInfos[baseNode].incidentCellsCount; incidentCellIndex++)
  {
    IndexType cellIndex = nodeTopologyInfos[baseNode].incidentCells[incidentCellIndex];
    for(IndexType incidentEdgeIndex = 0; incidentEdgeIndex < Space::EdgesPerCell; incidentEdgeIndex++)
    {
      IndexType edgeIndex = cellTopologyInfos[cellIndex].incidentEdges[incidentEdgeIndex];

      bool found = true;
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
      {
        if (nodeIndices[nodeNumber] != edges[edgeIndex].incidentNodes[nodeNumber]) found = false;
      }
      if (found)
      {
        return edgeIndex;
      }
    }
  }
  return IndexType(-1);//there is no such edge
}

template <typename Space>
typename Space::IndexType GeomMeshCommon<Space>::
  GetEdgeIndex(IndexType node0, IndexType node1) const
{
  IndexType nodeIndices[Space::NodesPerEdge];
  nodeIndices[0] = node0;
  nodeIndices[1] = node1;
  return GetEdgeIndex(nodeIndices);
}

template <>
bool GeomMeshCommon<Space3>::IsBoundaryCell(IndexType cellIndex) const
{
  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
  {
    IndexType faceIndex = cellTopologyInfos[cellIndex].incidentFaces[faceNumber];
    if (faceTopologyInfos[faceIndex].incidentCellsCount == 1) return true;
  }
  return false;
}

template <>
bool GeomMeshCommon<Space2>::IsBoundaryCell(IndexType cellIndex) const
{
  for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; ++edgeNumber)
  {
    IndexType edgeIndex = cellTopologyInfos[cellIndex].incidentEdges[edgeNumber];
    if (edgeTopologyInfos[edgeIndex].incidentCellsCount == 1) return true;
  }
  return false;
}

template <>
Space2::Scalar GeomMeshCommon<Space2>::GetAspectRatio(IndexType cellIndex) const
{
  Vector cellVertices[Space::NodesPerCell];
  GetCellVertices(cellIndex, cellVertices);

  Scalar area = Scalar(0.5) * fabs( (cellVertices[2] - cellVertices[0]) ^ (cellVertices[1] - cellVertices[0]) );

  Scalar perimeter = 0;
  Scalar maxEdge = 0;
  for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerCell; ++nodeNumber)
  {
    Vector edge = cellVertices[(nodeNumber + 1) % Space2::NodesPerCell] - cellVertices[nodeNumber];
    maxEdge = std::max(maxEdge, edge.Len());
    perimeter += edge.Len();
  }
  Scalar aspectRatio = maxEdge * perimeter * Scalar(0.25) / sqrt(3) / area;
  return aspectRatio;
}

template <>
Space3::Scalar GeomMeshCommon<Space3>::GetAspectRatio(IndexType cellIndex) const
{
  Vector cellVertices[Space::NodesPerCell];
  GetCellVertices(cellIndex, cellVertices);

  // http://docs.salome-platform.org/latest/gui/SMESH/aspect_ratio_3d_page.html

  Scalar maxEdge = 0;

  for (IndexType i0 = 0; i0 + 1 < Space::NodesPerCell; ++i0)
  {
    for (IndexType i1 = i0 + 1; i1 < Space::NodesPerCell; ++i1)
    {
      maxEdge = Max(maxEdge, (cellVertices[i1] - cellVertices[i0]).Len());
    }
  }

  Scalar alpha = fabs(MixedProduct(cellVertices[3] - cellVertices[0], cellVertices[2] - cellVertices[0], cellVertices[1] - cellVertices[0]));

  Scalar d = ((cellVertices[2] - cellVertices[0]) ^ (cellVertices[1] - cellVertices[0])).Len() +
    ((cellVertices[3] - cellVertices[0]) ^ (cellVertices[1] - cellVertices[0])).Len() +
    ((cellVertices[3] - cellVertices[0]) ^ (cellVertices[2] - cellVertices[0])).Len() +
    ((cellVertices[3] - cellVertices[1]) ^ (cellVertices[2] - cellVertices[1])).Len();

  return maxEdge / alpha * d / (2 * sqrt(Scalar(6.0)));
}

/*
template <>
bool GeomMeshCommon<Space2>::IsBoundaryNode(IndexType nodeIndex) const
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

template <>
bool GeomMeshCommon<Space3>::IsBoundaryNode(IndexType nodeIndex) const
{
  // TODO
  assert(0);
  return false;
}
*/

template <typename Space>
void GeomMeshCommon<Space>::BuildOriginalCellIndices(IndexType* originalIndices, IndexType cellsCount)
{
  originalCellIndices.resize(cellsCount);
  for (IndexType originalCellIndex = 0; originalCellIndex < cellsCount; ++originalCellIndex)
  {
    IndexType cellNodeIndices[Space::NodesPerCell];
    std::copy(originalIndices + originalCellIndex * Space::NodesPerCell, 
      originalIndices + originalCellIndex * Space::NodesPerCell + Space::NodesPerCell, cellNodeIndices);
    IndexType updatedCellIndex = GetCellIndex(cellNodeIndices);
    assert(updatedCellIndex != IndexType(-1));
    originalCellIndices[updatedCellIndex] = originalCellIndex;
  }
}

template <typename Space>
void GeomMeshCommon<Space>::BuildUpdatedCellIndices(IndexType* originalIndices, IndexType cellsCount)
{
  updatedCellIndices.resize(cellsCount);
  for (IndexType originalCellIndex = 0; originalCellIndex < cellsCount; ++originalCellIndex)
  {
    IndexType cellNodeIndices[Space::NodesPerCell];
    std::copy(originalIndices + originalCellIndex * Space::NodesPerCell, 
      originalIndices + originalCellIndex * Space::NodesPerCell + Space::NodesPerCell, cellNodeIndices);
    IndexType updatedCellIndex = GetCellIndex(cellNodeIndices);
    updatedCellIndices[originalCellIndex] = updatedCellIndex;
  }
}

template <typename Space>
void GeomMeshCommon<Space>::GetCellsPairOrientation(
  IndexType firstNodes[Space::NodesPerCell], IndexType secondNodes[Space::NodesPerCell],
  IndexType associatedPermutation[Space::NodesPerCell])
{
  for (IndexType i = 0; i < Space::NodesPerCell; ++i)
  {
    associatedPermutation[i] = IndexType(-1);
    for (IndexType j = 0; j < Space::NodesPerCell; ++j)
    {
      if (firstNodes[i] == secondNodes[j]) associatedPermutation[i] = j;
    }
    assert(associatedPermutation[i] != IndexType(-1));
  }
}

template <typename Space>
typename GeomMeshCommon<Space>::NodeGroupSize 
  GeomMeshCommon<Space>::GetNodeGroup(IndexType nodeIndex, std::vector<IndexType>& groupNodes) const
{
  assert(GetIncidentCellsCount(nodeIndex) > 0);
  IndexType incidentCellIndex = GetIncidentCellIndex(nodeIndex, 0);

  // TODO: replace with stack-memory-allocated structure
  std::stack<IndexType> s;
  std::set<IndexType> used;
  Vector node = nodes[nodeIndex].pos;

  s.push(incidentCellIndex);
  used.insert(incidentCellIndex);

  while (!s.empty())
  {
    IndexType currentCellIndex = s.top();
    s.pop();

    IndexType minNodeNumber = 0;
    for (IndexType nodeNumber = 1; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if ((node - nodes[cells[currentCellIndex].incidentNodes[nodeNumber]].pos).SquareLen() <
          (node - nodes[cells[currentCellIndex].incidentNodes[minNodeNumber]].pos).SquareLen())
      minNodeNumber = nodeNumber;
    }

    groupNodes.push_back(cells[currentCellIndex].incidentNodes[minNodeNumber]);

    Scalar minHeight = GetMinHeight(currentCellIndex);

    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
    {
      IndexType correspondingCellIndex = GetCorrespondingCellIndex(currentCellIndex, faceNumber);
      if (correspondingCellIndex != IndexType(-1))
      {
        bool isIncidentCellForNode = false;
        for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
        {
          IndexType currentNodeIndex = cells[correspondingCellIndex].incidentNodes[nodeNumber];
          if ((nodes[currentNodeIndex].pos - node).Len() < minHeight) isIncidentCellForNode = true;
        }

        if (isIncidentCellForNode && used.find(correspondingCellIndex) == used.end())
        {
          s.push(correspondingCellIndex);
          used.insert(correspondingCellIndex);
        }
      }
    }
  } 

  return groupNodes.size();
}

template <>
Space2::Scalar GeomMeshCommon<Space2>::GetVolume(IndexType cellIndex) const
{
  Vector cellVertices[Space::NodesPerCell];
  GetCellVertices(cellIndex, cellVertices);
  return Scalar(0.5) * fabs((cellVertices[2] - cellVertices[0]) ^ (cellVertices[1] - cellVertices[0]));
}

template <>
Space3::Scalar GeomMeshCommon<Space3>::GetVolume(IndexType cellIndex) const
{
  Vector cellVertices[Space::NodesPerCell];
  GetCellVertices(cellIndex, cellVertices);
  return fabs(MixedProduct(cellVertices[3] - cellVertices[0], cellVertices[2] - cellVertices[0], cellVertices[1] - cellVertices[0])) / Scalar(6.0);
}

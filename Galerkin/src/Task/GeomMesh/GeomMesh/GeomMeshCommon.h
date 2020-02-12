#pragma once

#include "GeomMeshBase.h"
#include "../../../Maths/Spaces.h"
#include "AdditionalCellInfo.h"
#include "../AABBTree.h"

#include <assert.h>
#include <vector>

template <typename Space>
struct GeomMeshCommon: public GeomMeshBase<Space>
{
  SPACE_TYPEDEFS
  typedef typename GeomMeshBase<Space>::CellTopologyInfo CellTopologyInfo;
  typedef IndexType NodeGroupSize;

  // geometry
  struct Node
  {
    SPACE_TYPEDEFS
    Node();
    Node(const Vector& pos);
    Vector pos;
  };

  struct Cell
  {
    SPACE_TYPEDEFS

    Cell();
    Cell(IndexType nodeIndices[Space::NodesPerCell]);

    bool operator==(const Cell& other) const;
    bool operator< (const Cell& other) const;
    bool operator!=(const Cell& other) const;

    Cell(const Cell& other);
    Cell& operator=(const Cell& other);

    IndexType incidentNodes[Space::NodesPerCell];

    friend std::ostream& operator<<(std::ostream& os, const Cell& cell)
    {
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
      {
        os << cell.incidentNodes[nodeNumber] << " ";
      }
      return os;
    }
    friend std::istream& operator>>(std::istream& os, Cell& cell)
    {
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
      {
        os >> cell.incidentNodes[nodeNumber];
      }
      return os;
    }
  private:
    void Copy(const Cell& other);
  };

  struct Edge
  {
    SPACE_TYPEDEFS
    IndexType incidentNodes[Space::NodesPerEdge];

    Edge();
    Edge(IndexType nodeIndices[Space::NodesPerEdge]);

    Edge(IndexType nodeIndex0, IndexType nodeIndex1);
    
    bool operator==(const Edge& other) const;
    bool operator< (const Edge& other) const;
    bool operator!=(const Edge& other) const;

    Edge(const Edge& other);
    Edge& operator=(const Edge& other);
  private:
    void Copy(const Edge& other);
  };

  Edge*       GetEdge(IndexType index);
  const Edge* GetEdge(IndexType index) const;

  Node*       GetNode(IndexType index);
  const Node* GetNode(IndexType index) const;

  Cell*       GetCell(IndexType index);
  const Cell* GetCell(IndexType index) const;

  IndexType GetNodesCount();
  IndexType GetEdgesCount();
  IndexType GetCellsCount();

  IndexType GetCellIndex(IndexType nodeIndices[Space::NodesPerCell]) const;
  IndexType GetEdgeIndex(IndexType node0, IndexType node1) const;
  IndexType GetEdgeIndex(IndexType nodeIndices[Space::NodesPerEdge]) const;

  void GetCellEdgeNodes(const IndexType* cellIncidentNodes, IndexType edgeNumber, IndexType* edgeNodes) const;
  void GetCellEdgeNodes(const IndexType cellIndex, IndexType edgeNumber, IndexType* edgeNodes) const;

  std::vector<Node> nodes;
  std::vector<Edge> edges;
  std::vector<Cell> cells;

  void LoadGeom(Vector* vertices, IndexType* indices, IndexType verticesCount, IndexType cellsCount);

  void GetCellVertices(IndexType cellIndex, Vector* points) const;
  void GetCellVertices(const IndexType* cellIncidentNodes, Vector* points) const;
  void GetFixedCellIndices(IndexType cellIndex, IndexType* incidentNodes) const;
  void GetCellIndices(IndexType cellIndex, IndexType* incidentNodes) const; // bad orientation is possible

  // mass center of all nodes
  Vector GetMassCenter() const;
  Vector GetMassCenter(IndexType cellIndex) const;

  Scalar GetMinHeight(IndexType cellIndex) const;
  Scalar GetMinHeight() const;

  Scalar GetAspectRatio(IndexType cellIndex) const;
  Scalar GetVolume(IndexType cellIndex) const;

  bool IsBoundaryCell(IndexType cellIndex) const;
  // bool IsBoundaryNode(IndexType nodeIndex) const;

  // topology
  struct NodeTopologyInfo
  {
    NodeTopologyInfo();
    IndexType  incidentCellsCount;
    IndexType* incidentCells;
  };

  struct EdgeTopologyInfo
  {
    EdgeTopologyInfo();
    IndexType* incidentCells;
    IndexType  incidentCellsCount;
  };

  std::vector<NodeTopologyInfo> nodeTopologyInfos;
  std::vector<EdgeTopologyInfo> edgeTopologyInfos;
  std::vector<CellTopologyInfo> cellTopologyInfos;

  IndexType GetIncidentCellsCount(IndexType nodeIndex) const;
  IndexType GetIncidentCellIndex(IndexType nodeIndex, IndexType cellNumber) const;
  NodeGroupSize GetNodeGroup(IndexType nodeIndex, std::vector<IndexType>& groupNodes) const;

  void BuildTopologyInfos();
  void BuildNodeTopologyInfos();
  void BuildEdgeTopologyInfos();
  void BuildCellTopologyInfos();

  std::vector<IndexType> nodeIncidentCellsPool;
  std::vector<IndexType> edgeIncidentCellsPool;
  std::vector<IndexType> cellIncidentEdgesPool;

  struct EdgeIndices
  {
    SPACE_TYPEDEFS
    
    EdgeIndices();
    EdgeIndices(IndexType nodeIndex0, IndexType nodeIndex1);

    bool operator< (const EdgeIndices& other) const;
    bool operator==(const EdgeIndices& other) const;
    bool operator!=(const EdgeIndices& other) const;


    IndexType nodeIndices[Space::NodesPerEdge];

    EdgeIndices(const EdgeIndices& other);
    EdgeIndices& operator=(const EdgeIndices& other);
    friend std::ostream& operator<<(std::ostream& os, const EdgeIndices& edgeIndices)
    {
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
      {
        os << edgeIndices.nodeIndices[nodeNumber] << " ";
      }
      return os;
    }
  private:
    void Copy(const EdgeIndices& other);
  };

  struct EdgePairIndices
  {
    EdgeIndices edges[2];
    bool operator==(const EdgePairIndices& other) const;
    bool IsEdgesCoincident() const;
    EdgePairIndices& operator=(const EdgePairIndices& other)
    {
      std::copy(other.edges, other.edges + 2, edges);
      return *this;
    }
    EdgePairIndices(const EdgePairIndices& other)
    {
      std::copy(other.edges, other.edges + 2, edges);
    }
    EdgePairIndices() {}
  };

  struct CellIndices
  {
    SPACE_TYPEDEFS
    IndexType nodeIndices[Space::NodesPerCell];
  };

  void GetCellsPairOrientation(IndexType firstNodes[Space::NodesPerCell], 
    IndexType secondNodes[Space::NodesPerCell],
    IndexType associatedPermutation[Space::NodesPerCell]);

  // additional cell info
  std::vector< AdditionalCellInfo<Space> > additionalCellInfos;
  IndexType GetCorrespondingCellIndex( IndexType cellIndex, IndexType faceNumber) const;
  IndexType GetCorrespondingFaceNumber(IndexType cellIndex, IndexType faceNumber) const;
  IndexType GetInteractionType(        IndexType cellIndex, IndexType faceNumber) const;

  std::vector<IndexType> originalCellIndices; // map from new cell indices to original
  std::vector<IndexType> updatedCellIndices;  // map from old cell indices to original 
  void BuildOriginalCellIndices(IndexType* originalIndices, IndexType cellsCount);
  void BuildUpdatedCellIndices(IndexType* originalIndices, IndexType cellsCount);

  // aabb tree
  virtual void BuildAABBTree();
  virtual void UpdateAABBTree();
  AABBTree<Space, IndexType> aabbTree;
  std::vector<IndexType> treeNodeCellIndices;
  AABB GetAABB() const;
  AABB GetCellAABB(IndexType cellIndex) const;
  void RemoveCellFromAABBTree(IndexType cellIndex);
  void AddToAABBTree(IndexType cellIndex);
  bool IsCellInAABBTree(IndexType cellIndex) const;
};

#include "GeomMeshCommon.inl"

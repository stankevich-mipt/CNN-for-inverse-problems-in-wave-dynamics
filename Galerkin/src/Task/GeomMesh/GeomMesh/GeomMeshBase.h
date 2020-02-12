#pragma once

#include "../../../Maths/Spaces.h"
#include "../TopologyReconstructor.h"

#include <assert.h>
#include <vector>

template <typename Space>
struct GeomMeshBase;

template <>
struct GeomMeshBase<Space2>
{
  virtual ~GeomMeshBase() = default;

  SPACE2_TYPEDEFS
  typedef Space2 Space;

  struct EdgeLocation
  {
    SPACE2_TYPEDEFS
    EdgeLocation();
    EdgeLocation(IndexType cellIndex, IndexType edgeNumber);

    bool IsNull() const;
    bool operator == (const EdgeLocation& other) const;
    bool operator != (const EdgeLocation& other) const;

    IndexType cellIndex;
    IndexType edgeNumber;
  };

  struct EdgeLocationPair
  {
    EdgeLocationPair();
    EdgeLocationPair(const EdgeLocation& firstEdge, const EdgeLocation& secondEdge);
    EdgeLocation edges[2];
  };

  struct CellTopologyInfo
  {
    IndexType incidentEdges[Space::EdgesPerCell];
  };
};

template <>
struct GeomMeshBase<Space3>
{
  SPACE3_TYPEDEFS
  typedef Space3 Space;

  GeomMeshBase(): submeshesCount(0), 
    submeshInfos(nullptr),
    faceBorderIndicesCount(0),
    edgeBorderIndicesCount(0),
    nodeBorderIndicesCount(0),
    faceContactIndicesCount(0),
    edgeContactIndicesCount(0),
    nodeContactIndicesCount(0),
    contactNodeGroupsCount(0),
    contactFaceGroupsCount(0),
    contactEdgeGroupsCount(0)
  {
  }

  struct CellTopologyInfo
  {
    IndexType incidentEdges[Space::EdgesPerCell];
    IndexType incidentFaces[Space::EdgesPerCell];
  };

  struct Face
  {
    SPACE3_TYPEDEFS
    IndexType  incidentNodes[Space3::NodesPerFace];
  };

  struct FaceTopologyInfo
  {
    SPACE3_TYPEDEFS
    IndexType* incidentCells;
    IndexType  incidentCellsCount;
  };

  struct FaceIndices
  {
    SPACE3_TYPEDEFS
    IndexType nodeIndices[Space3::NodesPerFace];
    bool operator==(const FaceIndices& other) const
    {
      for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
      {
        if (nodeIndices[nodeNumber] != other.nodeIndices[nodeNumber]) return false;
      }
      return true;
    }
    bool operator!=(const FaceIndices& other) const
    {
      return !this->operator==(other);
    }

    FaceIndices& operator=(const FaceIndices& other)
    {
      std::copy(other.nodeIndices, other.nodeIndices + Space::NodesPerFace, nodeIndices);
      return *this;
    }

    FaceIndices(const FaceIndices& other)
    {
      std::copy(other.nodeIndices, other.nodeIndices + Space::NodesPerFace, nodeIndices);
    }

    FaceIndices() {}

    friend std::ostream& operator<<(std::ostream& os, const FaceIndices& faceIndices)
    {
      for (Space::IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
      {
        os << faceIndices.nodeIndices[nodeNumber] << " ";
      }
      return os;
    }
  };

  struct BoundaryFace: public FaceIndices
  {
    BoundaryFace& operator=(const BoundaryFace& other)
    {
      std::copy(other.nodeIndices, other.nodeIndices + Space::NodesPerFace, nodeIndices);
      return *this;
    }
  };

  struct FacePairIndices
  {
    FaceIndices faces[2];
    
    FacePairIndices()
    {}

    FacePairIndices(const FacePairIndices& other)
    {
      Copy(other);
    }

    bool operator==(const FacePairIndices& other) const
    {
      return faces[0] == other.faces[0] && faces[1] == other.faces[1];
    }
    FacePairIndices& operator=(const FacePairIndices& other)
    {
      Copy(other);
      return *this;
    }

  private:
    void Copy(const FacePairIndices& other)
    {
      for (IndexType pairIndex = 0; pairIndex < 2; ++pairIndex)
      {
        faces[pairIndex] = other.faces[pairIndex];
      }
    }
  };

  struct NodePair
  {
    SPACE3_TYPEDEFS
    IndexType nodes[2];

    bool operator ==(const NodePair& p) const
    {
      return nodes[0] == p.nodes[0] && nodes[1] == p.nodes[1];
    }
    bool operator < (const NodePair& p) const
    {
      if (nodes[0] <  p.nodes[0])                          return true;
      if (nodes[0] == p.nodes[0] && nodes[1] < p.nodes[1]) return true;
      return false;
    }
  };

  struct FaceLocation
  {
    SPACE3_TYPEDEFS
    FaceLocation();
    FaceLocation(IndexType cellIndex, IndexType faceNumber);

    bool IsNull() const;
    bool operator == (const FaceLocation& other) const;
    bool operator != (const FaceLocation& other) const;

    IndexType cellIndex;
    IndexType faceNumber;
  };

  struct FaceLocationPair
  {
    FaceLocationPair();
    FaceLocationPair(const FaceLocation& firstFace, const FaceLocation& secondFace);
    FaceLocation faces[2];
  };

  struct EdgePair
  {
    SPACE3_TYPEDEFS

    IndexType edges[2];
    IndexType orientation;

    bool operator == (const EdgePair& p) const
    {
      return edges[0] == p.edges[0] && edges[1] == p.edges[1];
    }
    bool operator < (const EdgePair &p) const
    {
      if (edges[0] <  p.edges[0])                          return true;
      if (edges[0] == p.edges[0] && edges[1] < p.edges[1]) return true;
      return false;
    }
  };

  struct FacePair
  {
    SPACE3_TYPEDEFS
    IndexType faces[2];
    IndexType orientation;
  };

  struct ContactNode
  {
    SPACE3_TYPEDEFS
    bool operator < (const ContactNode& n) const
    {
      if (nodeIndex < n.nodeIndex) return true;
      return false;
    }
    IndexType nodeIndex;
  };

  struct ContactEdge
  {
    SPACE3_TYPEDEFS
    bool operator < (const ContactEdge& e) const
    {
      if (edgeIndex < e.edgeIndex) return true;
      return false;
    }
    IndexType incidentNodes[2];
    IndexType edgeIndex;
  };

  struct ContactFace
  {
    SPACE3_TYPEDEFS
    bool operator < (const ContactFace& f) const
    {
      if (faceIndex < f.faceIndex) return true;
      return false;
    }
    IndexType incidentNodes[3];
    IndexType faceIndex;
  };

  struct ContactNodePairInfo
  {
    SPACE3_TYPEDEFS
    Vector normal;
    IndexType typeIndex;
  };

  struct ContactEdgePairInfo
  {
    SPACE3_TYPEDEFS
    Vector normals[2];
    IndexType typeIndex;
  };

  struct ContactFacePairInfo
  {
    SPACE3_TYPEDEFS
    Vector normals[3];
    IndexType typeIndex;
  };

  std::vector<IndexType> nodeContactGroupIndices;
  std::vector<IndexType> edgeContactGroupIndices;
  std::vector<IndexType> faceContactGroupIndices;

  struct SubmeshInfo
  {
    IndexType firstNodeIndex;
    IndexType nodesCount;

    IndexType firstEdgeIndex;
    IndexType edgesCount;

    IndexType firstFaceIndex;
    IndexType facesCount;

    IndexType firstCellIndex;
    IndexType cellsCount;
  };

  IndexType submeshesCount;
  SubmeshInfo *submeshInfos;

  std::vector<Face> faces;
  std::vector<FaceTopologyInfo> faceTopologyInfos;

  std::vector<IndexType> faceIncidentCellsPool;

  std::vector<IndexType> faceBorderToGlobal;
  std::vector<IndexType> nodeGlobalToBorder;

  IndexType faceBorderIndicesCount;
  IndexType edgeBorderIndicesCount;
  IndexType nodeBorderIndicesCount;

  std::vector<Vector>    nodeBorderNormals;
  std::vector<FacePair>  faceContactToGlobal;
  std::vector<IndexType> nodeGlobalToContact;

  IndexType faceContactIndicesCount;
  IndexType edgeContactIndicesCount;
  IndexType nodeContactIndicesCount;

  struct ContactGroupInfo
  {
    IndexType groupInfoOffset;
    IndexType pairInfoOffset;
    IndexType size;
  };

protected:
  DuplicateRemover<IndexType, IndexType> edgeBorderToGlobal;
  DuplicateRemover<IndexType, IndexType> nodeBorderToGlobal;
  DuplicateRemover<EdgePair, IndexType>  edgeContactToGlobal;
  DuplicateRemover<NodePair, IndexType>  nodeContactToGlobal;

  struct ContactNodeGroupInfo: public ContactGroupInfo
  {
    IndexType borderInfoOffset;
    IndexType contactInfoOffset;
  };

  IndexType                     contactNodeGroupsCount;
  std::vector<ContactNode>      contactNodeGroupData;
  std::vector<ContactGroupInfo> contactNodeGroupInfos;

  IndexType                     contactFaceGroupsCount;
  std::vector<ContactFace>      contactFaceGroupData;
  std::vector<ContactGroupInfo> contactFaceGroupInfos;

  IndexType                     contactEdgeGroupsCount;
  std::vector<ContactEdge>      contactEdgeGroupData;
  std::vector<ContactGroupInfo> contactEdgeGroupInfos;

  std::vector<ContactNodePairInfo> contactNodePairInfos;
  std::vector<ContactEdgePairInfo> contactEdgePairInfos;
  std::vector<ContactFacePairInfo> contactFacePairInfos;
};

#include "GeomMeshBase.inl"



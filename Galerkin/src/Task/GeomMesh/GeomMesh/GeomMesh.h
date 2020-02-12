#pragma once

#include "../../../Maths/Spaces.h"
#include "GeomMeshCommon.h"

#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <algorithm>

template <typename Space>
struct GeomMesh;

template <>
struct GeomMesh<Space2>: public GeomMeshCommon<Space2>
{
  SPACE2_TYPEDEFS
  typedef Space2 Space;

  // geometry
  using GeomMeshCommon<Space>::GetCellIndex;
  using GeomMeshCommon<Space>::GetCorrespondingCellIndex;
  using GeomMeshCommon<Space>::GetCorrespondingFaceNumber;
  using GeomMeshCommon<Space>::GetInteractionType;
  
  IndexType GetCellIndex(IndexType node0, IndexType node1, IndexType node2) const;

  struct BoundaryEdge: public EdgeIndices
  {
    BoundaryEdge();
    BoundaryEdge(IndexType nodeIndex0, IndexType nodeIndex1);
  };

  void GetCellEdgeVertices(const IndexType cellIndex, IndexType edgeNumber, Vector* edgeVertices) const;
  Vector GetCellEdgeMiddle(const IndexType cellIndex, IndexType edgeNumber) const;

  void GetGhostCellVertices(IndexType cellIndex, IndexType boundaryEdgeNumber, Vector* ghostCellVertices) const;
  Vector GetEdgeExternalNormal(IndexType cellIndex, IndexType edgeNumber) const;
  Vector GetFaceExternalNormal(IndexType cellIndex, IndexType faceNumber) const; // it`s for convenience

  void GetCellFaceNodes(IndexType cellIndex, IndexType faceNumber, IndexType* faceNodes) const; // it`s for convenience

  bool             FindCommonEdge(IndexType srcCellIndex, IndexType dstCellIndex, IndexType& srcEdgeNumber, IndexType& dstEdgeNumber);
  IndexType        GetEdgeNumber(IndexType cellIndex, IndexType nodeIndex0, IndexType nodeIndex1) const;
  EdgeLocationPair BuildEdgeLocation(IndexType nodeIndex1, IndexType nodeIndex2);
  EdgeLocationPair BuildEdgeLocation(const EdgePairIndices& contactEdgePairIndices);

  void BuildAdditionalTopology(
    EdgePairIndices* contactEdges,  IndexType* contactEdgesCount,  IndexType contactTypesCount,
    BoundaryEdge*    boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount,
    IndexType *internalContactTypes = nullptr);

  EdgeLocationPair GetEdgeLocation(IndexType nodeIndex0, IndexType nodeIndex1) const;
  EdgeLocationPair GetEdgeLocation(EdgeIndices edgeIndices) const;

  bool IsBoundaryCell(IndexType cellIndex) const;
  bool IsBoundaryNode(IndexType nodeIndex) const;

  void FindNodeGroup(IndexType nodeIndex, IndexType* groupNodeIndices, IndexType& groupNodesCount);

private:
  void BuildCellAdditionalTopology(IndexType *internalContactTypes);
  void BuildContactsInfo(EdgePairIndices* contactEdges, IndexType* contactEdgesCount, IndexType contactTypesCount);
  void BuildBoundariesInfo(BoundaryEdge* boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount);
};

#include "GeomMesh2.inl"

template<>
class GeomMesh<Space3>: public GeomMeshCommon<Space3>
{
public:
  SPACE3_TYPEDEFS
  typedef Space3 Space;

  using GeomMeshCommon<Space>::GetCellIndex;
  using GeomMeshCommon<Space>::GetCorrespondingCellIndex;
  using GeomMeshCommon<Space>::GetCorrespondingFaceNumber;
  using GeomMeshCommon<Space>::GetInteractionType;

  void BuildAdditionalGroupsInfo(IndexType* submeshNodesCount, IndexType submeshesCount,
    FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
    BoundaryFace*    boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);

  void BuildAdditionalTopology(
    FacePairIndices* contactFaces,  IndexType* contactEdgesCount,  IndexType contactTypesCount,
    BoundaryFace*    boundaryFaces, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount,
    IndexType *internalContactTypes = nullptr);

  IndexType GetNodeSubmeshIndex(IndexType nodeIndex);
  IndexType GetEdgeSubmeshIndex(IndexType edgeIndex);
  IndexType GetFaceSubmeshIndex(IndexType faceIndex);
  IndexType GetCellSubmeshIndex(IndexType cellIndex);

  Vector GetFaceOuterNormal(IndexType  faceIndex,     Vector* vertexPositions);
  Vector GetFaceOuterNormal(IndexType* incidentNodes, Vector* vertexPositions);

  IndexType GetCellInNodeDirection(IndexType nodeIndex, Vector direction);
  IndexType GetCellInEdgeDirection(IndexType edgeIndex, Vector direction);
  IndexType GetCellInFaceDirection(IndexType faceIndex, Vector direction);

  IndexType GetCellIndex(IndexType node0, IndexType node1, IndexType node2, IndexType node3) const;

  IndexType GetFaceIndex(IndexType node0, IndexType node1, IndexType node2) const;
  IndexType GetFaceIndex(const IndexType nodeIndices[Space::NodesPerFace]) const;

  Face*       GetFace(IndexType index);
  const Face* GetFace(IndexType index) const;

  IndexType GetSubmeshNodesCount(IndexType submeshIndex);
  IndexType GetSubmeshEdgesCount(IndexType submeshIndex);

  IndexType GetFacesCount();
  IndexType GetSubmeshFacesCount(IndexType submeshIndex);

  IndexType GetSubmeshCellsCount(IndexType submeshIndex);
  IndexType GetSubmeshesCount();

  IndexType            GetContactNodeGroupsCount();
  IndexType            GetContactNodeGroupSize  (IndexType groupIndex);
  ContactNode          GetContactNodeInGroup    (IndexType groupIndex, IndexType sideIndex);
  IndexType            GetNodeSideInGroup       (IndexType groupIndex, IndexType nodeIndex);
  ContactNodePairInfo& GetContactNodePairInfo   (IndexType groupIndex, IndexType sideIndex0, IndexType sideIndex1);

  IndexType            GetContactEdgeGroupsCount();
  IndexType            GetContactEdgeGroupSize  (IndexType groupIndex);
  ContactEdge          GetContactEdgeInGroup    (IndexType groupIndex, IndexType sideIndex);
  IndexType            GetEdgeSideInGroup       (IndexType groupIndex, IndexType edgeIndex);
  ContactEdgePairInfo& GetContactEdgePairInfo   (IndexType groupIndex, IndexType sideIndex0, IndexType sideIndex1);

  IndexType            GetContactFaceGroupsCount();
  IndexType            GetContactFaceGroupSize  (IndexType groupIndex);
  ContactFace          GetContactFaceInGroup    (IndexType groupIndex, IndexType sideIndex);
  IndexType            GetFaceSideInGroup       (IndexType groupIndex, IndexType faceIndex);
  ContactFacePairInfo& GetContactFacePairInfo   (IndexType groupIndex, IndexType sideIndex0, IndexType sideIndex1);

  IndexType GetFacePairOrientation(IndexType srcCellNodes[3], IndexType dstCellNodes[3]);
  Scalar GetFaceSquare(Vector faceVertices[Space::NodesPerFace]);

  bool FindCommonFace(IndexType srcCellIndex, IndexType dstCellIndex, 
    IndexType &_srcFaceNumber, IndexType &_dstFaceNumber, IndexType &_orientation);

  void GetCellFaceNodes(const IndexType* cellIncidentNodes, IndexType faceNumber, IndexType* faceNodes) const;
  void GetCellFaceNodes(IndexType cellIndex, IndexType faceNumber, IndexType* faceNodes) const;
  void GetCellFaceVertices(IndexType cellIndex, IndexType faceNumber, Vector* faceVertices) const;

  void BuildBoundariesInfo(BoundaryFace* boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);
  IndexType GetLocalFaceNumber(IndexType cellIndex, IndexType* faceIndices);
  void BuildContactsInfo(FacePairIndices* contactFaces, IndexType* contactFacesCount, IndexType contactTypesCount);

  FaceLocationPair GetFaceLocation(IndexType nodeIndex0, IndexType nodeIndex1, IndexType nodeIndex2) const;
  FaceLocationPair GetFaceLocation(const IndexType nodeIndices[Space::NodesPerFace]) const;
  FaceLocationPair GetFaceLocation(const FacePairIndices& contactFacePairIndices) const;

  IndexType GetFaceNumber(IndexType cellIndex, IndexType nodeIndex0, IndexType nodeIndex1, IndexType nodeIndex2) const;
  IndexType GetFaceNumber(IndexType cellIndex, const IndexType faceNodeIndices[Space::NodesPerFace]) const;

  Vector GetFaceExternalNormal(IndexType cellIndex, IndexType faceNumber) const;
  Vector GetFaceExternalNormal(Vector* faceGlobalVertices) const;

  void GetGhostCellVertices(IndexType cellIndex, IndexType boundaryFaceNumber, Vector* ghostCellVertices) const;

private:
  void BuildSubmeshInfos(IndexType* submeshNodesCount, IndexType submeshesCount);

  void BuildFaceContactGroups(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                              FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);
  void BuildEdgeContactGroups(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                              FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);
  void BuildNodeContactGroups(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                              FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);

  void BuildContactNodePairInfos(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                                 FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);
  void BuildContactEdgePairInfos(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                                 FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);
  void BuildContactFacePairInfos(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                                 FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);
  void BuildSmoothContactNormals(FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                                 FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);
  void BuildHardContactNormals  (FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
                                 FaceIndices*     boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount);

  void ExpandContactPairInfos();

  void BuildCellAdditionalTopology(IndexType *internalContactTypes);
};

Space3::IndexType GetRelativeOrientationSide(Space3::IndexType* firstNodes, Space3::IndexType* secondNodes);
Space3::IndexType GetRelativeOrientationEdge(Space3::IndexType* firstNodes, Space3::IndexType* secondNodes);

bool PointIsInCell(Space3::Vector cellPoint0, 
                   Space3::Vector cellPoint1, 
                   Space3::Vector cellPoint2, 
                   Space3::Vector cellPoint3,
                   Space3::Vector testPoint);

#include "GeomMesh3.inl"


#pragma once

#include "../../../../Maths/Spaces.h"
#include "../../GeomMesh/GeomMesh.h"
#include <vector>
#include <numeric>

template <typename Space>
struct MeshIOBase;

template <>
struct MeshIOBase<Space2>
{
  SPACE2_TYPEDEFS

  typedef GeomMesh<Space2>::EdgePairIndices        EdgePairIndices;
  typedef GeomMesh<Space2>::BoundaryEdge           BoundaryEdge;

  // contacts description
  std::vector<IndexType> contactEdgesCount;
  std::vector<EdgePairIndices> contactEdges;

  // boundary description
  std::vector<IndexType> boundaryEdgesCount;
  std::vector<BoundaryEdge> boundaryEdges;

  bool operator==(const MeshIOBase<Space2>& other) const
  {
    return contactEdgesCount == other.contactEdgesCount &&
      contactEdges == other.contactEdges &&
      boundaryEdges == other.boundaryEdges &&
      boundaryEdgesCount == other.boundaryEdgesCount;
  }

  IndexType GetContactsCount() const
  {
    return std::accumulate(contactEdgesCount.begin(), contactEdgesCount.end(), 0);
  }

  IndexType GetBoundariesCount() const
  {
    return std::accumulate(boundaryEdgesCount.begin(), boundaryEdgesCount.end(), 0);
  }
};

template <>
struct MeshIOBase<Space3>
{
  SPACE3_TYPEDEFS

  typedef GeomMesh<Space3>::FacePairIndices FacePairIndices;
  typedef GeomMesh<Space3>::BoundaryFace    BoundaryFace;

  // contacts
  std::vector<IndexType>       contactFacesCount;
  std::vector<FacePairIndices> contactFaces;

  // boundary
  std::vector<IndexType>       boundaryFacesCount;
  std::vector<BoundaryFace>    boundaryFaces;

  bool operator==(const MeshIOBase<Space3>& other) const
  {
    bool countsEqual = (contactFacesCount == other.contactFacesCount && boundaryFacesCount == other.boundaryFacesCount);
    bool boundariesEqual = boundaryFaces == other.boundaryFaces;
    bool contactsEqual = true;

    return countsEqual && boundariesEqual && contactsEqual;
  }

  IndexType GetContactsCount() const
  {
    return std::accumulate(contactFacesCount.begin(), contactFacesCount.end(), 0);
  }

  IndexType GetBoundariesCount() const
  {
    return std::accumulate(boundaryFacesCount.begin(), boundaryFacesCount.end(), 0);
  }
};

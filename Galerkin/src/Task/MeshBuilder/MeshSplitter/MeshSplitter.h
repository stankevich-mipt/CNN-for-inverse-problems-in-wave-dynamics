#pragma once

#include "MeshSplitterCommon.h"

template <typename Space>
class MeshSplitter;

template <>
class MeshSplitter<Space2>: public MeshSplitterCommon<Space2>
{
public:
  SPACE2_TYPEDEFS
  typedef Space2                                                   Space;
  typedef GeomMesh<Space>                                          GeomMeshT;
  typedef GeomMeshT::BoundaryEdge                                  BoundaryEdge;
  typedef GeomMeshT::EdgePairIndices                               EdgePairIndices;
  typedef GeomMeshT::EdgeIndices                                   EdgeIndices;
  typedef GeomMeshT::EdgeLocationPair                              EdgeLocationPair;

  MeshSplitter(MeshIO<Space>* const mesh, GeomMesh<Space>* const geomMesh): 
    MeshSplitterCommon<Space>(mesh, geomMesh)
  {
  }

protected:
  void AddBoundariesToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains) override;
  void AddContactsToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains) override;
};

template <>
class MeshSplitter<Space3>: public MeshSplitterCommon<Space3>
{
public:
  SPACE3_TYPEDEFS
  typedef Space3                                                   Space;
  typedef GeomMesh<Space>                                          GeomMeshT;
  typedef GeomMeshT::BoundaryFace                                  BoundaryFace;
  typedef GeomMeshT::FacePairIndices                               FacePairIndices;
  typedef GeomMeshT::FaceIndices                                   FaceIndices;
  typedef GeomMeshT::FaceLocationPair                              FaceLocationPair;

  MeshSplitter(MeshIO<Space>* const mesh, GeomMesh<Space>* const geomMesh): 
    MeshSplitterCommon<Space>(mesh, geomMesh)
  {
  }

protected:
  void AddBoundariesToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains) override;
  void AddContactsToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains) override;
};

#include "MeshSplitter2.inl"
#include "MeshSplitter3.inl"

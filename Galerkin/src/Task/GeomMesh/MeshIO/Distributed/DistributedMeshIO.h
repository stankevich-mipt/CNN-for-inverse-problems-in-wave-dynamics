#pragma once

#include "../../../../Maths/Spaces.h"
#include "../Local/MeshIO.h"
#include "../../GeomMesh/GeomMesh.h"
#include "TransitionInfo.h"
#include <vector>
#include <fstream>

template <typename Space>
struct DistributedMeshIOCommon: public MeshIO<Space>
{
  SPACE_TYPEDEFS
  using MeshIO<Space>::GetCellsCount;
  using MeshIO<Space>::vertices;
  using MeshIO<Space>::indices; 
  using MeshIO<Space>::detectorsPositions;
  typedef typename GeomMesh<Space>::Cell Cell;
  typedef typename TransitionInfo<Space>::TransitionNode TransitionNode;

  IndexType index;
  IndexType domainsCount;

  std::vector< TransitionInfo<Space> > transitionInfos;

  DistributedMeshIOCommon(IndexType domainsCount);

  virtual void Save(const std::string& fileName, IO::FileType fileType = IO::Binary);
  virtual void Load(const std::string& fileName, IO::FileType fileType = IO::Binary);

  bool operator==(const DistributedMeshIOCommon& other) const;

  void SaveTransitionInfos(const std::string& vtkFileName);

private:
  void SaveTransitionInfo(std::fstream& file, IO::FileType fileType, const std::vector< TransitionInfo<Space> >&);
  void LoadTransitionInfo(std::fstream& file, IO::FileType fileType, std::vector< TransitionInfo<Space> >* const);
};

template <typename Space>
struct DistributedMeshIO;

template <>
struct DistributedMeshIO<Space2>: public DistributedMeshIOCommon<Space2>
{
  SPACE2_TYPEDEFS

  typedef GeomMesh<Space2>::EdgeIndices EdgeIndices;

  using MeshIO<Space2>::contactEdges;
  using MeshIO<Space2>::contactEdgesCount;
  using MeshIO<Space2>::contactTypesCount;
  using MeshIO<Space2>::boundaryEdges;
  using MeshIO<Space2>::boundaryEdgesCount;
  using MeshIO<Space2>::boundaryTypesCount;

  using DistributedMeshIOCommon<Space2>::transitionInfos;
  using DistributedMeshIOCommon<Space2>::SaveTransitionInfos;

  DistributedMeshIO(IndexType domainsCount): DistributedMeshIOCommon<Space2>(domainsCount)
  {}
};

template <>
struct DistributedMeshIO<Space3>: public DistributedMeshIOCommon<Space3>
{
  SPACE3_TYPEDEFS

  typedef GeomMesh<Space3>::FaceIndices FaceIndices;

  using MeshIO<Space3>::contactFaces;
  using MeshIO<Space3>::contactFacesCount;
  using MeshIO<Space3>::contactTypesCount;
  using MeshIO<Space3>::boundaryFaces;
  using MeshIO<Space3>::boundaryFacesCount;
  using MeshIO<Space3>::boundaryTypesCount;

  using DistributedMeshIOCommon<Space3>::transitionInfos;
  using DistributedMeshIOCommon<Space3>::SaveTransitionInfos;

  DistributedMeshIO(IndexType domainsCount): DistributedMeshIOCommon<Space3>(domainsCount)
  {}
};

#include "DistributedMeshIO.inl"

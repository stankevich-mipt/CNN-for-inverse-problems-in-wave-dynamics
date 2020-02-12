#pragma once

#include "../../../Maths/Spaces.h"
#include "../../GeomMesh/GeomMesh/AdditionalCellInfo.h"
#include "../../../Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO.h"
#include "../../GeomMesh/GeomMesh/GeomMesh.h"
#include "../../../Utils/Utils.h"

#include <vector>
#include <numeric>
#include <algorithm>

template <typename Space>
class MeshSplitterCommon
{
public:
  SPACE_TYPEDEFS
  typedef GeomMesh<Space>                                        GeomMeshT;
  typedef IndexType                                              LocalIndex;
  typedef IndexType                                              GlobalIndex;
  typedef typename AdditionalCellInfo<Space>::template AuxInfo<IndexType> CellInfo;
  typedef typename TransitionInfo<Space>::Cell Cell;
  typedef typename TransitionInfo<Space>::TransitionNode TransitionNode;

  struct Location
  {
    Location():
      domainIndex(IndexType(-1)),
      globalIndex(IndexType(-1))
    {}

    Location(IndexType domainIndex, IndexType objectIndex):
      domainIndex(domainIndex),
      globalIndex(objectIndex)
    {}
    IndexType domainIndex;

    union 
    {
      IndexType globalIndex;
      IndexType localIndex;
    };

    bool operator ==(const Location& other) const
    {
      return domainIndex == other.domainIndex && globalIndex == other.globalIndex;
    }

    bool operator < (const Location& other) const
    {
      if (domainIndex != other.domainIndex)
      {
        return domainIndex < other.domainIndex;
      } else
      {
        return globalIndex < other.globalIndex;
      }
    }
  };

  typedef std::map<Location, LocalIndex>  GlobalToLocalIndicesDictionary;

  MeshSplitterCommon(MeshIO<Space>* const mesh, GeomMesh<Space>* const geomMesh):
    mesh(mesh), geomMesh(geomMesh)
  {}

  virtual ~MeshSplitterCommon() = default;

  void Split(const std::vector<IndexType>& cellsDomainIds, 
    IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains, bool allowMovement);

  GlobalToLocalIndicesDictionary* GetCellsDictionary()
  {
    return &cellsDictionary;
  }

  // represents in which domain each detector locates
  std::vector<Location> detectorsLocations;

  struct DomainNodesIndices
  {
    std::vector<IndexType> globalIndices;
  };

  struct DomainCellsIndices
  {
    std::vector<IndexType> globalIndices;
  };

  struct DomainInfo
  {
    DomainInfo(IndexType domainsCount,
      IndexType contactTypesCount,
      IndexType boundaryTypesCount):
      contactsCount(contactTypesCount),
      boundariesCount(boundaryTypesCount)
    {}

    DomainNodesIndices nodesIndices;
    DomainCellsIndices cellsIndices;
    std::vector<IndexType> contactsCount;
    std::vector<IndexType> boundariesCount;
  };

  std::vector<DomainInfo> domainsInfos;
protected:
  // index of corresponding domain for each cell
  std::vector<IndexType> cellsDomainIds; 

  MeshIO<Space>* const mesh;
  GeomMeshT* geomMesh;

  // for each cell list of domains containing it
  std::vector<CellInfo> cellInfos;
  std::vector<IndexType> cellsCorrespondingDomainsPool;

  GlobalToLocalIndicesDictionary cellsDictionary;
  GlobalToLocalIndicesDictionary nodesDictionary;

  // computes mapping from global indices of vertices to reference ones
  void BuildNodesDictionary(IndexType domainsCount, bool allowMovement);
  void BuildCellCorrespondingDomainPool(IndexType domainsCount);

  // set correct sizes for data arrays
  void SetCorrectSizes(std::vector< DistributedMeshIO<Space> >* const domains);

  void AddCellToDomain(IndexType domainIndex, IndexType globalCellIndex, std::vector< DistributedMeshIO<Space> >* const domains);
  virtual void AddBoundariesToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains) = 0;
  virtual void AddContactsToDomains(IndexType domainsCount, std::vector< DistributedMeshIO<Space> >* const domains) = 0;
 
  void BuildTransitionCells(std::vector< DistributedMeshIO<Space> >* const domains);
  void BuildTransitionNodes(std::vector< DistributedMeshIO<Space> >* const domains);

  std::vector<IndexType> RenumerateCellsDomainIds(const std::vector<IndexType>& cellsDomainIds);

  struct TransitionCells
  {
    std::vector<IndexType> cellsIndices;
  };

  struct TransitionCellsDomainInfo
  {
    std::vector<TransitionCells> transitionCells;
  };

  std::vector<TransitionCellsDomainInfo> transitionCellDomainInfos;
};

#include "MeshSplitterCommon.inl"

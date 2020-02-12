#pragma once

#include "../Local/ElasticVolumeMesh.h"
#include "../../GeomMesh/MeshIO/Distributed/TransitionInfo.h"
#include "../../GeomMesh/MeshIO/Distributed/DistributedMeshIO.h"
#include "../../../DifferentialSolvers/DifferentialSolver.h"
#include <vector>
#include <map>

template<typename Space, typename FunctionSpace>
class DistributedElasticVolumeMeshCommon : public ElasticVolumeMesh<Space, FunctionSpace>
{
public:
  SPACE_TYPEDEFS

  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::CellSolution     CellSolution;
  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::Node             Node;
  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::Cell             Cell;
  typedef IndexType BytesHandled;

  const static int dimsCount      = ElasticVolumeMesh<Space, FunctionSpace>::ElasticSystemType::dimsCount;
  const static int functionsCount = FunctionSpace::functionsCount;

  DistributedElasticVolumeMeshCommon(DifferentialSolver<Scalar>* solver,
    Scalar tolerance,
    IndexType domainIndex,
    IndexType domainsCount,
    int hierarchyLevelsCount,
    bool sendNodesInfo) :
    ElasticVolumeMesh<Space, FunctionSpace>(solver, tolerance, hierarchyLevelsCount),
    domainIndex(domainIndex),
    domainsCount(domainsCount),
    syncDataSizes(domainsCount, 0),
    sendNodesInfo(sendNodesInfo)
  {
  }

  IndexType GetSyncDataSize(IndexType dstDomainIndex) const;

  using ElasticVolumeMesh<Space, FunctionSpace>::volumeMesh;
  using ElasticVolumeMesh<Space, FunctionSpace>::detectorsPositions;
  using ElasticVolumeMesh<Space, FunctionSpace>::MakeSnapshot;

  struct CellSyncData
  {
    CellSolution cellSolution;
    typename GeomMesh<Space>::Cell cell;
    bool destroyed;
  };

  struct NodeSyncData
  {
    Vector pos;
    IndexType nodeIndex;
  };

  void RebuildTimeHierarchyLevels(IndexType globalStepIndex, bool allowCollisions);
  bool IsSyncDataEmpty(IndexType dstDomainIndex) const;
  void SetTransitionInfo(const std::vector< TransitionInfo<Space> >& sendingInfo);

  void BuildCellsSyncData(IndexType dstDomainIndex, char* const cellsData);
  void UpdateCellsData(IndexType cellsCount, CellSyncData* const cells);

  void BuildNodesSyncData(IndexType dstDomainIndex, char* const nodesData);
  void UpdateNodesData(IndexType nodesCount, const NodeSyncData* const nodes);

  BytesHandled UpdateDomainData(const char* const data);
  BytesHandled BuildSyncData(IndexType dstDomainIndex, char* const data);
  virtual void ComputeSyncDataSizes();

protected:
  IndexType domainIndex;
  IndexType domainsCount;
  std::vector< TransitionInfo<Space> > transitionInfos;
  std::vector<IndexType> syncDataSizes;
  bool sendNodesInfo;

  typedef std::map<IndexType, IndexType> NodesDictionary;
  std::vector< NodesDictionary > nodesDictionaries;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include "DistributedElasticVolumeMeshCommon.inl"

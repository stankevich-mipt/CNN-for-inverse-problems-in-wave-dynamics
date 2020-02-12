#pragma once

#include "DistributedElasticVolumeMeshCommon.h"

template<typename Space, typename FunctionSpace>
class DistributedElasticVolumeMesh;

template<typename FunctionSpace>
class DistributedElasticVolumeMesh<Space2, FunctionSpace>: public DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>
{
public:
  SPACE2_TYPEDEFS
  typedef Space2 Space;

  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::EdgeLocationPair EdgeLocationPair;
  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::EdgeLocation     EdgeLocation;
  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::EdgeIndices      EdgeIndices;
  typedef typename DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::Node Node;
  typedef typename DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::Cell Cell;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::detectorsPositions;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainIndex;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::syncDataSizes;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::sendNodesInfo;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::MakeSnapshot;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::transitionInfos;

  DistributedElasticVolumeMesh(DifferentialSolver<Scalar>* solver, 
    Scalar tolerance,
    IndexType domainIndex,
    IndexType domainsCount,
    int hierarchyLevelsCount,
    bool sendNodesInfo = false,
    bool sendEdgesInfo = false):
  DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>(solver, 
    tolerance,
    domainIndex,
    domainsCount,
    hierarchyLevelsCount,
    sendNodesInfo),
    sendEdgesInfo(sendEdgesInfo)
  {
  }

  struct EdgeSyncData
  {
    typename GeomMesh<Space>::Edge edge;
    IndexType interactionType;
  };

  void UpdateDomainData(const char* const data);
  void BuildSyncData(IndexType dstDomainIndex, char* const data);

private:
  bool sendEdgesInfo;

  void BuildEdgesSyncData(IndexType dstDomainIndex, char* const data);
  void UpdateEdgesData(IndexType edgesCount, const EdgeSyncData* const edgesData);

  void ComputeSyncDataSizes() override;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#include "DistributedElasticVolumeMesh2.inl"

template <typename FunctionSpace>
class DistributedElasticVolumeMesh<Space3, FunctionSpace>: public DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>
{
public:
  SPACE3_TYPEDEFS
  typedef Space3 Space;

  using ElasticVolumeMesh<Space, FunctionSpace>::FindDestructions;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::detectorsPositions;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainIndex;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::syncDataSizes;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::sendNodesInfo;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::MakeSnapshot;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::transitionInfos;

  struct FaceSyncData
  {
    typename GeomMesh<Space>::Face face;
    IndexType interactionType;
  };

  DistributedElasticVolumeMesh(DifferentialSolver<Scalar>* solver, 
    Scalar tolerance,
    IndexType domainIndex,
    IndexType domainsCount,
    int hierarchyLevelsCount,
    bool sendNodesInfo = false,
    bool sendFacesInfo = false):
  DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>(solver, 
    tolerance,
    domainIndex,
    domainsCount,
    hierarchyLevelsCount,
    sendNodesInfo),
    sendFacesInfo(sendFacesInfo)
  {
  }

  void UpdateDomainData(const char* const data);
  void BuildSyncData(IndexType dstDomainIndex, char* const data);

private:
  bool sendFacesInfo;

  void BuildFacesSyncData(IndexType dstDomainIndex, char* const data);
  void UpdateFacesData(IndexType facesCount, const FaceSyncData* const facesData);

  void ComputeSyncDataSizes() override;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include "DistributedElasticVolumeMesh3.inl"

#pragma once

#include <limits>
#include <algorithm>
#include "../../ElasticSystem/ElasticSystem.h"
#include "../../ElasticSystem/IniStates.h"
#include "../../../IO/Vtk/BasicVtkWriter.h"
#include "../../../Maths/Spaces.h"
#include "../../VolumeMethod/VolumeMesh/VolumeMesh.h"
#include "../../Task/SettingsParser/SolverSettingsParser.h"
#include "../../GeomMesh/NodeGroupManager.h"

template<typename Space, typename FunctionSpace>
struct ElasticVolumeMeshCommon: public DifferentialSystem<typename Space::Scalar>
{
  SPACE_TYPEDEFS

  typedef          ElasticSystem<Space>                                ElasticSystemType;
  typedef typename ElasticSystemType::ElasticSpace                     ElasticSpaceType;
  typedef typename ElasticSystemType::MediumParameters                 MediumParameters;
  typedef          VolumeMesh<Space, FunctionSpace, ElasticSystemType> VolumeMeshType;
  typedef          VolumeMeshCommon<Space, FunctionSpace, ElasticSystemType> VolumeMeshTypeCommon;
  typedef typename ElasticSpaceType::Elastic                           Elastic;
  typedef typename VolumeMeshType::Node                                Node;
  typedef typename VolumeMeshType::Cell                                Cell;
  typedef typename VolumeMeshType::CellSolution                        CellSolution;
  typedef          Space                                               SpaceT;

  const static int functionsCount = FunctionSpace::functionsCount;
  const static int dimsCount      = VolumeMeshType::dimsCount;

  ElasticVolumeMeshCommon(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  std::vector<Vector> detectorsPositions;
  std::vector<IndexType> detectorsCells;

  VolumeMeshType volumeMesh;
  DifferentialSolver<Scalar>* solver;
  Scalar tolerance;

  IndexType unfoldIterationsCount;
  Scalar minGridHeight;

  Scalar tensionErrorMult;
  Scalar velocityErrorMult;
  Scalar positionErrorMult;

  Scalar tensionDimensionlessMult;
  Scalar velocityDimensionlessMult;

  bool allowMovement;
  bool allowPlasticity;
  bool allowContinuousDestruction;
  bool allowDiscreteDestruction;

  IndexType updateCollisionInfoPeriod;

  void Initialize(bool allowMovement,
    Scalar collisionWidth,
    IndexType unfoldIterationsCount, Scalar minGridHeight,
    Scalar tensionErrorMult,
    Scalar velocityErrorMult,
    Scalar positionErrorMult,
    Scalar tensionDimensionlessMult,
    Scalar velocityDimensionlessMult);

  void LoadState(const std::vector< IniStateMaker<ElasticSpaceType>* >& stateMakers, Scalar mult);

  Elastic InterpolateElasticRef(IndexType cellIndex, Vector refCoords) const;
  Elastic InterpolateElasticRef(IndexType cellIndex, IndexType pointIndex, typename VolumeMeshTypeCommon::PointType pointType) const;
  Elastic InterpolateElastic(IndexType cellIndex, Vector globalPoint, bool halfStepSolution = false) const;

  Elastic GetAverageCellElastic(IndexType cellIndex) const;
  
  int GetDimentionsCount(const SolverState& solverState) const override;
  int GetMaxDimentionsCount() const override;

  void GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState) override;
  void GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState& solverState) override;
  void GetCurrCoords(Scalar& time, Scalar* currCoords) const override;

  void SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState) override;
  void SetCurrCoords(Scalar time, const Scalar* newCoords, const Scalar* oldCoords, const SolverState& solverState) override;
  void SetCurrCoords(Scalar time, const Scalar* oldCoords) override;
  
  void SetDetectors(const std::vector<Vector>& globalDetectorsPositions);
  void GetDetectorsData(IndexType detectorIndex, Scalar* data);

  Scalar GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState&, const Scalar* mults) override;

  void FindDestructions(std::vector<bool>* isCellBroken);
  virtual void UnfoldMesh(Scalar minHeight, IndexType iterationsCount) = 0;
  virtual void HandleContinuousDestruction() = 0;

  void HandleDamping();
  void SetDamping(Scalar damping);
  Scalar GetDamping() const;

  void DestroyCell (IndexType cellIndex, Scalar powderShearMult);
  ElasticSystemType* GetSystem();
  void MoveSceneToSnapshotRegion();

  /* 
    sqrt(2/3 * strain_{ij} : strain_{ij})
  */
  Scalar GetDeformRate(IndexType cellIndex, Vector refCoords) const;

  // stress correction due to plasticity in compliance with  Prandtl-Ruess model
  void HandlePlasticity(Scalar dt);

  typename SolverSettings<Space>::Erosion erosion;
  typename SolverSettings<Space>::DynamicContactBox dynamicContactBox;

  void   GetCellEnergy(IndexType cellIndex, Scalar& kineticEnergy, Scalar& potentialEnergy) const;
  Vector GetTotalImpulse() const;
  Scalar GetTotalEnergy(Scalar& kineticEnergy, Scalar& potentialEnergy) const;
  Scalar GetTotalMass() const;

  void HandleMaterialErosion();

  virtual Scalar GetTimeStepPrediction() override;

  virtual void DestroyFace(IndexType cellIndex, IndexType faceNumber, IndexType dynamicBoundaryType) = 0;

  // total work of the plasticity for each cell
  std::vector<Scalar> plasticDeforms;

  std::vector<bool> isCellBroken;

  bool ProcessPlasticity(const Scalar k, Elastic& elastic, bool updateElastic = false);

  void MakeRhoCorrection(Scalar minTimeStep = 2e-5);

  void SetGlobalStepIndex(IndexType globalStepIndex);

protected:
  void ComputeElasticMults()
  {
    ComputeElasticMults(Overload<Space>());
  }

  IndexType globalStepIndex;

  void UpdateAABBTree(IndexType globalStepIndex);

  NodeGroupManager<Space, VolumeMeshType> nodeGroupManager;

  // for error computation
  Scalar elasticMults[dimsCount];
  Scalar damping;
  Scalar minHeightInMesh;

  struct PrecomputedBasisFunctions
  {
    Polynomial<Scalar, IndexType, Space::Dimension> basisFunctions[functionsCount];
    Polynomial<Scalar, IndexType, Space::Dimension> basisDerivativeFunctions[functionsCount][Space::Dimension];
  } precomputedBasisFunctions;

private:
  void ComputeElasticMults(Overload<Space2>)
  {
    for (int valueIndex = 0; valueIndex < volumeMesh.dimsCount; ++valueIndex)
    {
      Scalar mult = Scalar(1.0);
      switch (valueIndex)
      {
        case 0: case 1: case 2:
        {
          mult = tensionErrorMult;
        } break;
        case 3: case 4:
        {
          mult = velocityErrorMult;
        } break;
      }
      elasticMults[valueIndex] = mult;
    }
  }

  void ComputeElasticMults(Overload<Space3>)
  {
    for (int valueIndex = 0; valueIndex < volumeMesh.dimsCount; ++valueIndex)
    {
      Scalar mult = Scalar(1.0);
      switch (valueIndex)
      {
        case 0: case 1: case 2: case 3: case 4: case 5:
        {
          mult = tensionErrorMult;
        } break;
        case 6: case 7: case 8:
        {
          mult = velocityErrorMult;
        } break;
      }
      elasticMults[valueIndex] = mult;
    }
  }
protected:
  // for test only
  bool CheckNodeGroupInfo(IndexType nodeIndex)
  {
    IndexType nodeGroupSize = nodeGroupManager.GetGroupSize(nodeIndex);

    if (nodeGroupSize == 0) return true;

    std::vector<IndexType> pool;
    pool.reserve(60);

    int test[99] = { 0 };
    IndexType size = volumeMesh.GetNodeGroup(nodeIndex, pool);
    
    if (size != nodeGroupSize)
    {
      std::cout << nodeIndex << "\n";
      return false;
    }

    const IndexType* nodeGroupP = nodeGroupManager.GetGroup(nodeIndex);
    std::copy(nodeGroupP, nodeGroupP + size, test);

    std::sort(pool.begin(), pool.end());
    std::sort(test, test + size);

    for (int i = 0; i < size; ++i)
    {
      if (test[i] != pool[i])
      {
        return false;
      }
    }
    return true;
  }

  bool CheckFace(IndexType cellIndex, IndexType faceNumber)
  {
    IndexType faceIndices[Space::NodesPerFace];
    volumeMesh.GetCellFaceNodes(cellIndex, faceNumber, faceIndices);

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
    {
      if (!CheckNodeGroupInfo(faceIndices[nodeNumber])) return false;
    }
    return true;
  }

  bool CheckCell(IndexType cellIndex)
  {
    IndexType cellIndices[Space::NodesPerCell];
    volumeMesh.GetFixedCellIndices(cellIndex, cellIndices);

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (!CheckNodeGroupInfo(cellIndices[nodeNumber])) return false;
    }
    return true;
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename Space, typename FunctionSpace>
struct ElasticVolumeMesh;

template<typename FunctionSpace>
struct ElasticVolumeMesh<Space2, FunctionSpace>: public ElasticVolumeMeshCommon<Space2, FunctionSpace>
{
  SPACE2_TYPEDEFS
  typedef          Space2                                                  Space;
  typedef          ElasticSystem<Space>                                    ElasticSystemType;
  using ElasticSpaceType = ElasticSystemType::ElasticSpace;
  using Elastic = ElasticSpaceType::Elastic;
  typedef          VolumeMesh<Space, FunctionSpace, ElasticSystemType>     VolumeMeshType;
  typedef typename VolumeMeshType::EdgePairIndices                         EdgePairIndices;
  typedef typename VolumeMeshType::BoundaryEdge                            BoundaryEdge;
  typedef typename VolumeMeshType::EdgeLocation                            EdgeLocation;
  typedef typename VolumeMeshType::EdgeLocationPair                        EdgeLocationPair;
  typedef typename VolumeMeshType::EdgeIndices                             EdgeIndices;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::tensionDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::velocityDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetAverageCellElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetSystem;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCell;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::erosion;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::minHeightInMesh;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::HandleMaterialErosion;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowContinuousDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowDiscreteDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::FindDestructions;

  ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  virtual void MakeSnapshot(Elastic* destData,
    IndexType width, IndexType height, 
    Vector boxPoint1, Vector boxPoint2, bool halfStepSolution);

  virtual void MakeSnapshot(Elastic* destData,
    const Vector& origin, const Vector& spacing, 
    const Vector& boxPoint1, const Vector& boxPoint2,
    bool halfStepSolution);

  Elastic GetAverageEdgeElastic (IndexType cellIndex, IndexType edgeNumber) const;

  void HandleContinuousDestruction() override;

  void DestroyFace(IndexType cellIndex, IndexType edgeNumber, IndexType dynamicBoundaryType) override
  {
    IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, edgeNumber);
    IndexType correspondingEdgeNumber = volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber;

    volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType = dynamicBoundaryType;

    if (correspondingCellIndex != IndexType(-1) && correspondingEdgeNumber != IndexType(-1))
    {
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringEdges[correspondingEdgeNumber].correspondingCellIndex = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringEdges[correspondingEdgeNumber].correspondingEdgeNumber = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringEdges[correspondingEdgeNumber].interactionType = dynamicBoundaryType;
    }
    this->nodeGroupManager.UpdateFace(cellIndex, edgeNumber);
  }

private:
  void UnfoldMesh(Scalar minHeight, IndexType iterationsCount) override;
};

template <typename FunctionSpace>
struct ElasticVolumeMesh<Space3, FunctionSpace>: public ElasticVolumeMeshCommon<Space3, FunctionSpace>
{
  SPACE3_TYPEDEFS
  typedef          Space3                                       Space;
  typedef          ElasticSystem<Space3>                        ElasticSystemType;
  using ElasticSpaceType = ElasticSystemType::ElasticSpace;
  using MediumParameters = ElasticSystemType::MediumParameters;

  typedef          VolumeMesh<Space3, FunctionSpace, ElasticSystemType>    VolumeMeshType;
  typedef typename VolumeMeshType::FacePairIndices                         FacePairIndices;
  typedef typename VolumeMeshType::BoundaryFace                            BoundaryFace;
  typedef typename VolumeMeshType::Node                                    Node;
  typedef typename VolumeMeshType::Cell                                    Cell;
  typedef typename VolumeMeshType::CellSolution                            CellSolution;
  using Elastic = ElasticSpaceType::Elastic;
  typedef typename VolumeMeshType::FaceLocation                            FaceLocation;
  typedef typename VolumeMeshType::FaceLocationPair                        FaceLocationPair;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::tensionDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::velocityDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetAverageCellElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetSystem;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCell;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::HandlePlasticity;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::HandleMaterialErosion;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowContinuousDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowDiscreteDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::FindDestructions;

  ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  void UnfoldMesh(Scalar minHeight, IndexType iterationsCount) override;

  virtual void MakeSnapshot(Elastic* destData,
    const Vector& origin, const Vector& spacing, 
    const Vector& boxPoint1, const Vector& boxPoint2,
    bool halfStepSolution);

  void HandleContinuousDestruction() override;

  void DestroyFace(IndexType cellIndex, IndexType faceNumber, IndexType dynamicBoundaryType) override
  {
    IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);
    IndexType correspondingFaceNumber = volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber;

    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].orientation = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType = dynamicBoundaryType;

    if (correspondingCellIndex != IndexType(-1) && correspondingFaceNumber != IndexType(-1))
    {
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].correspondingCellIndex = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].correspondingFaceNumber = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].orientation = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].interactionType = dynamicBoundaryType;
    }
    this->nodeGroupManager.UpdateFace(cellIndex, faceNumber);
  }

};

template<typename MeshType, typename StateMaker>
struct FunctionGetter
{
  typedef typename MeshType::SpaceT Space;
  SPACE_TYPEDEFS

  FunctionGetter(MeshType* mesh, StateMaker* stateMaker, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->stateMaker = stateMaker;
    this->cellIndex = cellIndex;
    mesh->volumeMesh.GetCellVertices(cellIndex, cellVertices);
  }

  void operator()(const Vector& point, Scalar* values)
  {
    Vector globalPoint = mesh->volumeMesh.RefToGlobalVolumeCoords(point, cellVertices);

    Scalar lambda = mesh->volumeMesh.cellMediumParameters[cellIndex].lambda;
    Scalar mju    = mesh->volumeMesh.cellMediumParameters[cellIndex].mju;
    Scalar invRho = mesh->volumeMesh.cellMediumParameters[cellIndex].invRho;

    typename MeshType::Elastic elastic = stateMaker->GetValue(globalPoint, lambda, mju, invRho);
    for (IndexType valueIndex = 0; valueIndex < mesh->dimsCount; ++valueIndex)
    {
      values[valueIndex] = elastic.values[valueIndex];
    }
  }

  StateMaker* stateMaker;
  MeshType*   mesh;
  IndexType   cellIndex;
  Vector      cellVertices[Space::NodesPerCell];
};

template <typename MeshType, typename Corrector>
void ApplyCorrector(MeshType* const mesh)
{
  using Scalar = typename MeshType::Scalar;
  using IndexType = typename MeshType::IndexType;
  using Vector = typename MeshType::Vector;

  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(mesh->volumeMesh.cells.size()); ++cellIndex)
  {
    Corrector corrector(mesh, cellIndex);
    Scalar cellValues[mesh->functionsCount * mesh->dimsCount];
    
    mesh->volumeMesh.functionSpace->template DecomposePrecomputed< Corrector, MeshType::dimsCount >(corrector, cellValues);
    //mesh->volumeMesh.functionSpace->template Decompose< Corrector, MeshType::dimsCount >(corrector, cellValues);


    for (IndexType functionIndex = 0; functionIndex < mesh->functionsCount; ++functionIndex)
    {
      for (IndexType valueIndex = 0; valueIndex < mesh->dimsCount; ++valueIndex)
      {
        mesh->volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
          cellValues[valueIndex * mesh->functionsCount + functionIndex];
      }
    }
  }
}

#include "ElasticVolumeMeshCommon.inl"
#include "ElasticVolumeMesh2.inl"
#include "ElasticVolumeMesh3.inl"

#pragma once

#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <sstream>
#include "assert.h"
#include <omp.h>
#include <errno.h>
#include <ctime>

#include "../../Network/NetworkInterface.h"
#include "../../Utils/Utils.h"

#include "../VolumeMethod/Spaces/Decomposer.h"
#include "../VolumeMethod/Spaces/PolynomialPrecomputer.h"
#include "../VolumeMethod/Spaces/DubinerPrecomputer.h"

#include "../VolumeMethod/Spaces/PolynomialSpace.h"
#include "../VolumeMethod/Spaces/FourierSpace2.h"
#include "../VolumeMethod/Spaces/LagrangeSpace.h"

#include "../ElasticVolumeMesh/Distributed/DistributedElasticVolumeMesh.h"

#include "../../IO/Vtk/BasicVtkWriter.h"
#include "../../IO/Vtk/SnapshotVtkWriter.h"
#include "../../IO/Vtk/MeshVtkWriter.h"
#include "../../IO/Vtk/ContactVtkWriter.h"
#include "../../IO/Vtk/CellInfoVtkWriter.h"

#include "../../DifferentialSolvers/SolversFactory.h"

#include "../ElasticSystem/SourceFunctors.h"
#include "../ElasticSystem/VectorFunctors.h"
#include "SettingsParser/SettingsParser.h"
#include "../GeomMesh/MeshIO/Distributed/DistributedMeshIO.h"
#include "../../Maths/Spaces.h"


// #define PROFILING
// #define SINGLE_THREAD
// #define WRITE_ENERGY_AND_IMPULSE

template<typename Space, unsigned int order>
class Task: public NotifyListener, public ReceiveListener
{
public:
  SPACE_TYPEDEFS

  typedef          ElasticSystem<Space>                       ElasticSystemType;
  typedef typename ElasticSystemType::ElasticSpace            ElasticSpace;
  typedef typename ElasticSpace::MediumParametersType         MediumParameters;
  typedef typename ElasticSpace::Elastic                      Elastic;
  const static int dimsCount = ElasticSystemType::dimsCount;

  // typedef PolynomialPrecomputer<Space, QuadratureDecomposer<Space, LagrangeSpace<Space, order> > > FunctionSpace;
  typedef PolynomialPrecomputer<Space, QuadratureDecomposer<Space, PolynomialSpace<Space, order> > > FunctionSpace;
 // typedef PolynomialPrecomputer<Space, NodalDecomposer <Space, LagrangeSpace<Space, order> > > FunctionSpace;
  // typedef DubinerPrecomputer<Space, QuadratureDecomposer<Space, DubinerSpace<Space, order> > > FunctionSpace;

  Task();
  virtual ~Task();
  void Run();

protected:
  void BuildContactDescriptions();
  void BuildBoundaryDescriptions();

  void LoadMeshes();
  void SaveVtkSnapshot(Scalar currTime, IndexType snapshotIndex, IndexType stepIndex, bool halfStepSolution = false);

  void AnalyzeSnapshotData(Scalar currTime, IndexType snapshotIndex, IndexType stepIndex, Overload<Space2>);
  void AnalyzeSnapshotData(Scalar currTime, IndexType snapshotIndex, IndexType stepIndex, Overload<Space3>);
protected:
  void SaveDetectorsData(Scalar currTime, IndexType stepIndex);
  void SaveDetectorsDataToFile();

  void BuildSourceTerms();
  void BuildPointSources();
  void AllocateSnapshotData();
  void BuildIniStateMakers();
  void LoadInitialState();
  void FindDestructions(IndexType domainNumber);
  MediumParameters MakeElasticMediumParams(typename MeshSettings<Space>::MediumParamsSection::MediumParams params);

private:
  void SaveProfilingData(const std::string& info, double begin, double end); 
  void ComputePacketsToReceiveCount();
 
  void SetTransitionInfo(IndexType domainNumber, DistributedMeshIO<Space>* mesh)
  {
    distributedElasticMeshes[domainNumber]->SetTransitionInfo(mesh->transitionInfos);
  }

  void LoadGeom(IndexType domainNumber);

  void LoadGeom(Overload<Space2>, IndexType domainNumber)
  {
    distributedElasticMeshes[domainNumber]->volumeMesh.LoadGeom(
    meshes[domainNumber]->vertices.data(), meshes[domainNumber]->indices.data(), 
    meshes[domainNumber]->vertices.size(), meshes[domainNumber]->GetCellsCount(),
    meshes[domainNumber]->contactEdges.data(), 
    meshes[domainNumber]->contactEdgesCount.data(), 
    meshes[domainNumber]->contactTypesCount,
    meshes[domainNumber]->boundaryEdges.data(),
    meshes[domainNumber]->boundaryEdgesCount.data(),
    meshes[domainNumber]->boundaryTypesCount,
    cellMediumParams.data(),
    internalContactTypes.data());
  }
  void LoadGeom(Overload<Space3>, IndexType domainNumber)
  {
    distributedElasticMeshes[domainNumber]->volumeMesh.LoadGeom(
      meshes[domainNumber]->vertices.data(), meshes[domainNumber]->indices.data(),
      meshes[domainNumber]->vertices.size(), meshes[domainNumber]->GetCellsCount(),
      meshes[domainNumber]->contactFaces.data(),
      meshes[domainNumber]->contactFacesCount.data(),
      meshes[domainNumber]->contactTypesCount,
      meshes[domainNumber]->boundaryFaces.data(),
      meshes[domainNumber]->boundaryFacesCount.data(),
      meshes[domainNumber]->boundaryTypesCount,
      cellMediumParams.data(),
      internalContactTypes.data());
  }

  IndexType FindSnapshotSize()
  {
    IndexType maxSnapshotSize = 0;
    for (IndexType snapshotIndex = 0; snapshotIndex < settings.snapshots.size(); ++snapshotIndex)
    {
      IntAABB iaabb;
      snapshotsSizes[snapshotIndex] = ::ComputeSnapshotAABBSize<Space>(settings.snapshots[snapshotIndex].data.origin, 
        settings.snapshots[snapshotIndex].data.spacing, 
        AABB(settings.snapshots[snapshotIndex].data.boxPoint1, settings.snapshots[snapshotIndex].data.boxPoint2));
      maxSnapshotSize = std::max(maxSnapshotSize, snapshotsSizes[snapshotIndex].GetVolume());
    }
    return maxSnapshotSize;
  }

  VectorFunctor<Space>* CreateBoundaryFunctor(std::vector<typename MeshSettings<Space>::BoundarySection::VectorFunctor> &functors)
  {
    struct CombinedFunctor: public VectorFunctor<Space>
    {
      virtual Vector operator ()(const Vector& point, const Vector& norm, Scalar time) const override
      {
        Vector res = Vector::zero();
        for (size_t functorIndex = 0; functorIndex < functors.size(); functorIndex++)
        {
          res += functors[functorIndex]->operator()(point, norm, time);
        }
        return res;
      }
      void Add(VectorFunctor<Space> *newbie)
      {
        functors.push_back(newbie);
      }
      void SetCurrentVelocity(const Vector& v) override
      {
        for (IndexType functorIndex = 0; functorIndex < functors.size(); ++functorIndex)
        {
          functors[functorIndex]->SetCurrentVelocity(v);
        }
      }

      bool IsEmpty() const
      {
        return functors.empty();
      }
      std::vector<VectorFunctor<Space> *> functors;
    };
    CombinedFunctor *combinedFunctor = new CombinedFunctor();
    typedef typename MeshSettings<Space>::BoundarySection BoundarySection;
    for (size_t functorIndex = 0; functorIndex < functors.size(); functorIndex++)
    {
      typename MeshSettings<Space>::BoundarySection::VectorFunctor &func = functors[functorIndex];
      switch (func.type)
      {
        case BoundarySection::VectorFunctor::Radial:
        {
          typename BoundarySection::RadialFunctorInfo radialInfo = settings.mesh.boundarySection.radialFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new RadialVectorFunctor<Space>(radialInfo.center, radialInfo.intensity));
        }break;
        case BoundarySection::VectorFunctor::Box:
        {
          typename BoundarySection::BoxFunctorInfo boxInfo = settings.mesh.boundarySection.boxFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new BoxVectorFunctor<Space>(AABB(boxInfo.boxPoint1, boxInfo.boxPoint2), boxInfo.value, boxInfo.aabbVelocity));
        }break;
        case BoundarySection::VectorFunctor::Const:
        {
          typename BoundarySection::ConstFunctorInfo constInfo = settings.mesh.boundarySection.constFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new ConstVectorFunctor<Space>(constInfo.value));
        }break;
        case BoundarySection::VectorFunctor::Wave:
        {
          typename BoundarySection::WaveFunctorInfo waveInfo = settings.mesh.boundarySection.waveFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new WaveVectorFunctor<Space>(
            waveInfo.waveVector, waveInfo.waveLength, waveInfo.initialPhase, waveInfo.speed, waveInfo.maxTime, waveInfo.shear));
        }break;
        case BoundarySection::VectorFunctor::Rotating:
        {
          typename BoundarySection::RotatingFunctorInfo rotatingInfo = settings.mesh.boundarySection.rotatingFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new RotatingVectorFunctor<Space>(
            rotatingInfo.pos, rotatingInfo.angularVelocity, rotatingInfo.linearVelocity, rotatingInfo.linearTime));
        }break;
        case BoundarySection::VectorFunctor::Rieker:
        {
          typename BoundarySection::RiekerFunctorInfo riekerInfo = settings.mesh.boundarySection.riekerFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new RiekerVectorFunctor<Space>(riekerInfo.waveVector, riekerInfo.peakFrequency));
        }break;
        case BoundarySection::VectorFunctor::Tabulated:
        {
          typename BoundarySection::TabulatedFunctorInfo tabulatedInfo = settings.mesh.boundarySection.tabulatedFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new TabulatedVectorFunctor<Space>(tabulatedInfo.direction, tabulatedInfo.fileName));
        }break;
        case BoundarySection::VectorFunctor::HydraulicPressure:
        {
          typename BoundarySection::HydraulicPressureFunctorInfo hydraulicPressureInfo =
            settings.mesh.boundarySection.hydraulicPressureFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new HydraulicPressureFunctor<Space>(hydraulicPressureInfo.fluidRho,
            hydraulicPressureInfo.g, hydraulicPressureInfo.fluidSurfacePoint));
        } break;
        case BoundarySection::VectorFunctor::HydrodynamicResistance:
        {
          typename BoundarySection::HydrodynamicResistanceFunctorInfo hydrodynamicResistanceInfo =
            settings.mesh.boundarySection.hydrodynamicResistanceFunctorInfos[func.infoIndex];
          combinedFunctor->Add(new HydrodynamicResistanceFunctor<Space>(
            AABB(hydrodynamicResistanceInfo.boxPoint1, hydrodynamicResistanceInfo.boxPoint2), hydrodynamicResistanceInfo.mediumRho,
            hydrodynamicResistanceInfo.flowVelocity, hydrodynamicResistanceInfo.alpha, hydrodynamicResistanceInfo.beta));
        } break;
      }
    }
    return combinedFunctor->IsEmpty() ? 0 : combinedFunctor;
  }

  static void PrintSolverState(const SolverState& solverState, Scalar currTime, Scalar timeStep);

  void UpdateMeshData(char* data);
  void OnNotify() override;
  void OnDataReceive(void *data, int, int) override;
  void SynchronizeMeshes();

  void      LoadNodesSchedule();
  IndexType GetDomainsCount() const;
  IndexType GetNodeId() const;
  IndexType GetCurrentNodeDomainsCount() const;
  static void SetThreadsCount();

  void WriteEnergyAndImpulseDeviation(const std::vector<Scalar>& initialEnergies, 
    const std::vector<Vector>& initialImpulses, Scalar currTime);

protected:
  Settings<Space> settings;
  SolverState     solverState;
  std::vector< DifferentialSolver<Scalar>* > solvers;

  typedef DistributedElasticVolumeMesh<Space, FunctionSpace> DistributedElasticMesh;

  std::vector<DistributedElasticMesh* > distributedElasticMeshes;

  std::vector<typename ElasticSystem<Space>::MediumParameters> cellMediumParams;
  std::vector<IndexType> internalContactTypes;

  // structure for loading mesh from file
  std::vector<DistributedMeshIO<Space>* > meshes;
  // data for each edge
  std::vector< std::vector<bool> > isContactBroken;
  // data for each cell
  std::vector< IniStateMaker<ElasticSpace>* > stateMakers;

  struct Sample
  {
    Scalar time;
    Elastic elastic;
  };

  typedef std::vector< std::vector<Sample> > DomainDetectorsData;
  std::vector< DomainDetectorsData > detectorsData;

  // parameter of distribution
  // lists of domains which will be computed on the i-th computation node
  std::vector<typename ScheduleSettings<Space>::NodeSchedule> nodesSchedule;
  // computation node index for each domain
  std::vector<IndexType> domainsLocation;
  // local index of all current node domains (domainIndex -> domainNumber)
  std::vector<IndexType> domainsNumbers;

  // network
  NetworkInterface* network;

  // profiling
  std::fstream profilingDataFile;

  struct TransitionData
  {
    char* sendingData;
  };
  std::vector< std::vector<TransitionData> > syncData;

  IndexType         notificationsToReceive;
  IndexType         receivedPackets;
  IndexType         packetsToReceive;

  struct LogicSumComparator
  {
    bool operator()(bool oldValue, bool newValue) const
    {
      return oldValue && newValue;
    }
  } logicSumComparator;

  struct SumComparator
  {
    IndexType operator()(IndexType oldValue, IndexType newValue) const
    {
      return oldValue + newValue;
    }
  } sumComparator;

  struct MinValueComparator
  {
    Scalar operator()(Scalar oldValue, Scalar newValue)
    {
      return std::min<Scalar>(oldValue, newValue);
    }
  } minValueComparator;

  // saving
  typename ElasticSpace::Elastic* snapshotData;
  std::vector<IndexVector> snapshotsSizes;

  std::fstream snapshotCollectionFile;
  
  SnapshotVtkWriter<ElasticSpace>          snapshotWriter;

  typedef typename AdditionalCellInfo<Space>::template AuxInfo<char> CharCellInfo;
  MeshVtkWriter< Space, CharCellInfo >      meshWriter;

  ContactVtkWriter<Space, FunctionSpace> contactsWriter;
  CellInfoVtkWriter<Space, FunctionSpace> cellInfosWriter;
};

template<typename Space, unsigned int order>
void Task<Space, order>::Run()
{
  network->init(nullptr, nullptr);
  double taskBegin = MPI_Wtime();
  SetThreadsCount();

  printf("Node %d network initialized\n", static_cast<int>(network->getID()));
  network->setReceiveListener(this);
  settings.Parse("task.xml");

  printf("Starting task with config %s in %dd space\n", settings.settingsFileName.c_str(), settings.configDimsCount);
  assert(settings.configDimsCount == Space::Dimension);
  LoadNodesSchedule();

  char profilingFileName[100];
  sprintf(profilingFileName, "out/profiling[%d].txt", int(GetNodeId()));
  profilingDataFile.open(profilingFileName, std::fstream::out);

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    solvers.push_back(SolversFactory<Scalar, IndexType>::Build(settings.solver.integrator));
  }

  distributedElasticMeshes.resize(GetCurrentNodeDomainsCount());
  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    distributedElasticMeshes[domainNumber] = new DistributedElasticMesh(
      solvers[domainNumber],
      settings.solver.tolerance,
      nodesSchedule[GetNodeId()].domainsIndices[domainNumber],
      GetDomainsCount(), settings.solver.hierarchyLevelsCount, 
      settings.solver.allowMovement,
      settings.solver.allowDiscreteDestruction || settings.solver.allowContinuousDestruction);
    distributedElasticMeshes[domainNumber]->allowPlasticity  = settings.solver.allowPlasticity;
    distributedElasticMeshes[domainNumber]->allowContinuousDestruction = settings.solver.allowContinuousDestruction;
    distributedElasticMeshes[domainNumber]->allowDiscreteDestruction   = settings.solver.allowDiscreteDestruction;

    distributedElasticMeshes[domainNumber]->updateCollisionInfoPeriod = settings.solver.updateCollisionInfoPeriod;

    distributedElasticMeshes[domainNumber]->erosion = settings.solver.erosion;
    distributedElasticMeshes[domainNumber]->dynamicContactBox = settings.solver.dynamicContactBox;
  }
  detectorsData.resize(GetCurrentNodeDomainsCount());

  LoadMeshes();
  BuildSourceTerms();
  BuildPointSources();
  BuildIniStateMakers();
  LoadInitialState();

  snapshotsSizes.resize(settings.snapshots.size());
  AllocateSnapshotData();

  {
    char filename[1024];
    sprintf(filename, "out/SnapshotCollection[%d].pvd" , static_cast<int>(GetNodeId()));
    snapshotCollectionFile.open(filename, fstream::out);
  }

  if (snapshotCollectionFile.fail())
  {
    std::cerr << "Can`t open SnapshotCollection.pvd";
    throw;
  } else
  {
    snapshotWriter.WriteSnapsotCollectionHeader(snapshotCollectionFile);
  }

  Scalar timeStep = settings.solver.maxTimeStep * settings.solver.maxScale;
  bool forceStep  = false;
  std::vector<Scalar> lastSnapshotTimes(settings.snapshots.size(), 0.0);
  Scalar currTime = Scalar(0.0);

  // all meshes have the same phases count
  solverState.phasesCount = distributedElasticMeshes[0]->solver->GetPhasesCount();
  solverState.hierarchyLevelsCount = settings.solver.hierarchyLevelsCount;

  int lastInitialGlobalStep = -1;
  Scalar lastInitialCurrTime = Scalar(-1);
  const IndexType MaxHierarchyLevel = 1 << (settings.solver.hierarchyLevelsCount - 1);
  bool phaseAdvanced = true;

  printf("Entering main cycle\n");
  PrintSolverState(solverState, currTime, timeStep);

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    distributedElasticMeshes[domainNumber]->RebuildTimeHierarchyLevels(solverState.globalStepIndex, settings.solver.allowMovement);
    distributedElasticMeshes[domainNumber]->solver->InitStep(timeStep, settings.solver.tolerance, true);
    distributedElasticMeshes[domainNumber]->solver->InitStep(solverState);
    distributedElasticMeshes[domainNumber]->SetDamping(exp(-settings.solver.damping * timeStep));
  }

  double startStep = MPI_Wtime();

  std::vector<Scalar> initialEnergies;
  std::vector<Vector> initialImpulses;
  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    Scalar kineticEnergy, potentialEnergy;
    initialEnergies.push_back(distributedElasticMeshes[domainNumber]->GetTotalEnergy(kineticEnergy, potentialEnergy));
    initialImpulses.push_back(distributedElasticMeshes[domainNumber]->GetTotalImpulse());
  }

  while(currTime < settings.task.destinationTime || settings.task.destinationTime < 0)
  {
    if (phaseAdvanced && solverState.AllHierarchyLevelsCompleted())
    {
      if (lastInitialCurrTime < currTime && lastInitialGlobalStep < solverState.globalStepIndex)
      {
        double savingBegin = MPI_Wtime();
        SaveDetectorsData(currTime, solverState.globalStepIndex / MaxHierarchyLevel);        
        bool snapshot = false;
        for (IndexType snapshotIndex = 0; snapshotIndex < settings.snapshots.size(); ++snapshotIndex)
        {
          IndexType currentStepIndex = (solverState.globalStepIndex / MaxHierarchyLevel);
          if ((currTime >= lastSnapshotTimes[snapshotIndex] + settings.snapshots[snapshotIndex].timePeriod && 
               currTime >= settings.snapshots[snapshotIndex].startTime &&
               currTime <= settings.snapshots[snapshotIndex].finishTime) || 
              (currentStepIndex % settings.snapshots[snapshotIndex].framePeriod == 0 && 
               currentStepIndex >= settings.snapshots[snapshotIndex].startFrame && 
               currentStepIndex <= settings.snapshots[snapshotIndex].finishFrame))
          {
            lastSnapshotTimes[snapshotIndex] = currTime;
            SaveVtkSnapshot(currTime, snapshotIndex, solverState.globalStepIndex);

            //For task-specific snapshot analyzing and debugging purposes
            //AnalyzeSnapshotData(currTime, snapshotIndex, solverState.globalStepIndex, Overload<Space>());

            snapshot = true;
          }
        }
        if (snapshot) 
        {
          SaveDetectorsDataToFile();
          WriteEnergyAndImpulseDeviation(initialEnergies, initialImpulses, currTime);
        }
        lastInitialCurrTime = currTime;
        lastInitialGlobalStep = solverState.globalStepIndex;
        double savingEnd = MPI_Wtime();
        SaveProfilingData(" Saving snapshot", savingBegin, savingEnd);
      }
    }
     
    bool allPhasesCompleted;
    SolverState nextState = solverState.GetNext(&allPhasesCompleted);
    if (!solverState.IsUseful())
    {
      solverState = nextState;
      continue;
    }

    printf(".");
    phaseAdvanced = true;
    for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
    {
      double beginPhase = MPI_Wtime();
      phaseAdvanced = phaseAdvanced && distributedElasticMeshes[domainNumber]->solver->AdvancePhase(solverState);
      double endPhase = MPI_Wtime();
      SaveProfilingData("   Advance phase", beginPhase, endPhase);
    }

    if (!allPhasesCompleted || !phaseAdvanced)
    {
      SynchronizeMeshes();
      if (phaseAdvanced)
      {
        solverState = nextState;
      }
    } else
    {
      double startAllPhasesCompleted = MPI_Wtime();

      bool globalStepSuccessful = true;
      printf("\n");
      for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
      {
        Scalar localError = distributedElasticMeshes[domainNumber]->solver->GetLastStepError();
        printf("node %d; domain %d; error %f\n", static_cast<int>(GetNodeId()), static_cast<int>(nodesSchedule[GetNodeId()].domainsIndices[domainNumber]), localError);

        bool stepSuccessful = forceStep || (localError < Scalar(2.0));

        if(!stepSuccessful)
        {
          printf("domain %d at node %d error is too large, recomputing global step\n",
            (int)nodesSchedule[GetNodeId()].domainsIndices[domainNumber], (int)GetNodeId());
        }
        // notify everyone about stepSuccessful received from others
        globalStepSuccessful = network->template Negotiate<bool, LogicSumComparator>(stepSuccessful, logicSumComparator);
      }

      if (globalStepSuccessful)
      {
        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
        {
          distributedElasticMeshes[domainNumber]->SetGlobalStepIndex(solverState.globalStepIndex);
          distributedElasticMeshes[domainNumber]->solver->AdvanceStep(solverState);
        }
        currTime = distributedElasticMeshes[0]->solver->GetCurrTime();

        double beginDestruction = MPI_Wtime();
        if (solverState.IsPreInitial())
        {
          for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
          {
            if (settings.solver.allowContinuousDestruction || settings.solver.allowDiscreteDestruction) 
              FindDestructions(domainNumber);
            if (settings.mesh.moveMassCenter)     distributedElasticMeshes[domainNumber]->MoveSceneToSnapshotRegion();
            if (settings.solver.allowPlasticity)  distributedElasticMeshes[domainNumber]->HandlePlasticity(
              distributedElasticMeshes[domainNumber]->solver->GetCurrStep());
            if (settings.solver.damping > 0)      distributedElasticMeshes[domainNumber]->HandleDamping();
          }
        }
        double endDestruction = MPI_Wtime();
        SaveProfilingData("   Destruction", beginDestruction, endDestruction);

        SynchronizeMeshes();
        solverState = nextState;
      } else
      {
        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
        {
          distributedElasticMeshes[domainNumber]->solver->RevertStep(lastInitialCurrTime);
        }
        solverState.SetInitialState(lastInitialGlobalStep);
        currTime = lastInitialCurrTime;
      }

      if (!globalStepSuccessful || solverState.AllHierarchyLevelsCompleted())
      {
        Scalar globalDesiredStep = std::numeric_limits<Scalar>::max() * Scalar(0.5);

        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
        {
          Scalar localDesiredStep = distributedElasticMeshes[domainNumber]->solver->GetTimeStepPrediction();
          localDesiredStep = minValueComparator(globalDesiredStep, localDesiredStep);
          // notify everyone about localDesiredStep received from others
          globalDesiredStep = network->template Negotiate<Scalar, MinValueComparator>(localDesiredStep, minValueComparator);
        }

        timeStep = globalDesiredStep;
        forceStep = false;

        if (timeStep > settings.solver.maxTimeStep) timeStep = settings.solver.maxTimeStep;
        if (timeStep < settings.solver.maxTimeStep * settings.solver.maxScale)
        {
          timeStep = settings.solver.maxTimeStep * settings.solver.maxScale;
          forceStep = true;
        }

        for (IndexType snapshotIndex = 0; snapshotIndex < settings.snapshots.size(); ++snapshotIndex)
        {
          if(currTime + timeStep * MaxHierarchyLevel > 
              lastSnapshotTimes[snapshotIndex] + settings.snapshots[snapshotIndex].timePeriod + std::numeric_limits<Scalar>::epsilon())
          {
            Scalar desiredStep = fabs(lastSnapshotTimes[snapshotIndex] + settings.snapshots[snapshotIndex].timePeriod - currTime) / MaxHierarchyLevel;
            desiredStep = std::max(desiredStep, settings.solver.maxTimeStep * settings.solver.maxScale / MaxHierarchyLevel);
            timeStep = std::min(timeStep, desiredStep);
            printf("Time step truncated by snapshot time\n");
          }
        }

        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
        {
          if (settings.solver.allowContinuousDestruction || settings.solver.allowDiscreteDestruction || 
            settings.solver.allowMovement || settings.solver.allowPlasticity)
          {
            distributedElasticMeshes[domainNumber]->RebuildTimeHierarchyLevels(solverState.globalStepIndex, settings.solver.allowMovement);
          }
          distributedElasticMeshes[domainNumber]->solver->InitStep(timeStep, settings.solver.tolerance, globalStepSuccessful);
          distributedElasticMeshes[domainNumber]->SetDamping(exp(-settings.solver.damping * timeStep));
        }

      }

      for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
      {
        distributedElasticMeshes[domainNumber]->solver->InitStep(solverState);
      }

      PrintSolverState(solverState, currTime, timeStep);
      double endStep = MPI_Wtime();
      SaveProfilingData(" AllPhasesCompleted", startAllPhasesCompleted, endStep);
      SaveProfilingData("Step", startStep, endStep);

      startStep = MPI_Wtime();
    }
  }
  SaveDetectorsDataToFile();
  snapshotWriter.WriteSnapsotCollectionLastTags(snapshotCollectionFile);
  snapshotCollectionFile.close();
  network->finalize();
  double taskEnd = MPI_Wtime();
  std::cout << "Total time elapsed: " << taskEnd - taskBegin;
}

template<typename Space, unsigned int order>
void Task<Space, order>::BuildContactDescriptions()
{
  typedef typename MeshSettings<Space>::ContactSection ContactSection;
  typedef typename ContactSection::Contact Contact;

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    ElasticSystem<Space>* system = distributedElasticMeshes[domainNumber]->GetSystem();
    ContactSection contactSection = settings.mesh.contactSection;

    for(IndexType contactIndex = 0; contactIndex < contactSection.contacts.size(); ++contactIndex)
    {
      IndexType interactionType = contactSection.contacts[contactIndex].interactionType;
      switch(contactSection.contacts[contactIndex].type)
      {
        case Contact::Glue:
        {
          IndexType infoIndex = contactSection.contacts[contactIndex].infoIndex;
          if(infoIndex != IndexType(-1))
          {
            typename ContactSection::GlueInfo glueInfo = contactSection.glueInfos[infoIndex];
            system->SetGlueContact(interactionType, glueInfo.maxShearStress / settings.solver.tensionDimensionlessMult, 
                                                    glueInfo.maxLongitudinalStress / settings.solver.tensionDimensionlessMult, 
                                                    glueInfo.dynamicBoundaryInteractionType);
          }else
          {
            system->SetGlueContact(interactionType);
          }
        }break;
        case Contact::Glide:
        {
          system->SetGlideContact(interactionType);
        }break;
        case Contact::Boundary:
        {
          IndexType infoIndex = contactSection.contacts[contactIndex].infoIndex;
          if(infoIndex != IndexType(-1))
          {
            typename ContactSection::BoundaryInfo boundaryInfo = contactSection.boundaryInfos[infoIndex];
            //system->SetBoundaryContact(interactionType, boundaryInfo.boundaryInteractionType);
          }
        }break;
        case Contact::Friction:
        {
          IndexType infoIndex = contactSection.contacts[contactIndex].infoIndex;
          if(infoIndex != IndexType(-1))
          {
            typename ContactSection::FrictionInfo frictionInfo = contactSection.frictionInfos[infoIndex];
            system->SetFrictionContact(interactionType, frictionInfo.frictionCoeff, frictionInfo.dynamicBoundaryInteractionType);
          }else
          {
            system->SetFrictionContact(interactionType);
          }
        }break;
      }
    }
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::BuildBoundaryDescriptions()
{
  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    ElasticSystem<Space> *system = distributedElasticMeshes[domainNumber]->GetSystem();
    typename MeshSettings<Space>::BoundarySection boundarySection = settings.mesh.boundarySection;
    typedef typename MeshSettings<Space>::BoundarySection BoundarySection;
    typedef typename MeshSettings<Space>::BoundarySection::Boundary Boundary;
    for(IndexType boundaryIndex = 0; boundaryIndex < boundarySection.boundaries.size(); ++boundaryIndex)
    {
      VectorFunctor<Space> *functor = nullptr;

      IndexType interactionType = boundarySection.boundaries[boundaryIndex].interactionType;
      Scalar    reflectionCoeff = boundarySection.boundaries[boundaryIndex].reflectionCoeff;

      switch(boundarySection.boundaries[boundaryIndex].type)
      {
        case Boundary::Free:
        {
          IndexType infoIndex = boundarySection.boundaries[boundaryIndex].infoIndex;
          if(infoIndex != IndexType(-1))
          {
            IndexType dynamicContactInteractionType = boundarySection.freeInfos[infoIndex].dynamicContactInteractionType;
            functor = CreateBoundaryFunctor(boundarySection.freeInfos[infoIndex].forceFunctors);
            system->SetFreeBoundary(interactionType, settings.solver.tensionDimensionlessMult, reflectionCoeff, functor, dynamicContactInteractionType);
          }else
          {
            system->SetFreeBoundary(interactionType, settings.solver.tensionDimensionlessMult, reflectionCoeff);
          }
        }break;
        case Boundary::Fixed:
        {
          IndexType infoIndex = boundarySection.boundaries[boundaryIndex].infoIndex;
          if(infoIndex != IndexType(-1))
          {
            functor = CreateBoundaryFunctor(boundarySection.fixedInfos[infoIndex].velocityFunctors);
            system->SetFixedBoundary(interactionType, settings.solver.velocityDimensionlessMult, reflectionCoeff, functor);
          }else
          {
            system->SetFixedBoundary(interactionType, settings.solver.velocityDimensionlessMult, reflectionCoeff);
          }
        }break;
        case Boundary::Absorb:
        {
          system->SetAbsorbBoundary(interactionType);
        }break;
        case Boundary::Symmetry:
        {
          system->SetSymmetryBoundary(interactionType);
        }break;
        case Boundary::AntiSymmetry:
        {
          system->SetAntiSymmetryBoundary(interactionType);
        }break;
        default:
          assert(0);
        break;
      }
    }
  }
}

template<typename Space, unsigned int order>
typename Task<Space, order>::MediumParameters Task<Space, order>::MakeElasticMediumParams(
  typename MeshSettings<Space>::MediumParamsSection::MediumParams params)
{
  typename ElasticSystem<Space>::MediumParameters res;
  res.lambda  = params.lambda;
  res.mju     = params.mju;
  res.invRho  = Scalar(1.0) / params.rho;

  res.initialRho = params.rho;

  res.plasticity.k0 = params.k;
  res.plasticity.a = params.alpha;
  res.plasticity.brittle = params.brittle;
  res.fixed = params.fixed;
  res.plasticity.maxPlasticDeform = params.maxPlasticDeform;
  res.plasticity.powderShearMult = params.powderShearMult;

  // grad(flowVelocity) must be equiled 0. it`s standart divergence-free condition
  res.flowVelocity = params.flowVelocity;
  res.MakeDimensionless(settings.solver.tensionDimensionlessMult, settings.solver.velocityDimensionlessMult);
  return res;
}

template<typename Space, unsigned int order>
void Task<Space, order>::FindDestructions(IndexType domainNumber)
{
  distributedElasticMeshes[domainNumber]->FindDestructions(&distributedElasticMeshes[domainNumber]->isCellBroken);
}

template<typename Space, unsigned int order>
void Task<Space, order>::SaveVtkSnapshot(Scalar currTime, 
  IndexType snapshotIndex, IndexType stepIndex, bool halfStepSolution)
{
  Vector boxPoint1 = settings.snapshots[snapshotIndex].data.boxPoint1;
  Vector boxPoint2 = settings.snapshots[snapshotIndex].data.boxPoint2;
  Vector origin = settings.snapshots[snapshotIndex].data.origin;
  Vector spacing = settings.snapshots[snapshotIndex].data.spacing;

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    // snapshot
    // AABB2 snapshotArea = distributedElasticMeshes[domainNumber]->volumeMesh.GetComputationArea();
    distributedElasticMeshes[domainNumber]->MakeSnapshot(snapshotData, 
      origin / settings.solver.velocityDimensionlessMult, 
      spacing / settings.solver.velocityDimensionlessMult, 
      boxPoint1 / settings.solver.velocityDimensionlessMult, 
      boxPoint2 / settings.solver.velocityDimensionlessMult, halfStepSolution);


    IndexType domainIndex = nodesSchedule[GetNodeId()].domainsIndices[domainNumber];

    char stepString[256];
    char domainString[256];
    char timeString[256];
    sprintf(stepString, "%.6d", static_cast<int>(stepIndex));
    sprintf(domainString, "%.3d", static_cast<int>(domainIndex));
    sprintf(timeString, "%.5f", currTime);

    if(settings.snapshots[snapshotIndex].data.used)
    {
      std::string dataFilename = settings.snapshots[snapshotIndex].data.filename;
      ReplaceSubstring(dataFilename, "<step>",    std::string(stepString));
      ReplaceSubstring(dataFilename, "<domain>",  std::string(domainString));
      ReplaceSubstring(dataFilename, "<time>",    std::string(timeString));

      printf("saving snapshot, %s\n\n", dataFilename.c_str());

      assert(dataFilename.find(".vti") != std::string::npos);
      snapshotWriter.Write(dataFilename.c_str(),
        snapshotData,
        origin, spacing,
        AABB(boxPoint1, boxPoint2),
        AABB(boxPoint1, boxPoint2),
        settings.snapshots[snapshotIndex].data.writeVelocity,
        settings.snapshots[snapshotIndex].data.writeTension,
        settings.snapshots[snapshotIndex].data.writeWaves);
      snapshotCollectionFile << "<DataSet timestep=\"" << currTime << "\" group=\"\" part=\"0\" file=\"" << dataFilename << "\"/>\n";
      snapshotCollectionFile.flush();
    }

    if(settings.snapshots[snapshotIndex].mesh.used)
    {
      std::string meshFilename = settings.snapshots[snapshotIndex].mesh.filename;
      ReplaceSubstring(meshFilename, "<step>",    std::string(stepString));
      ReplaceSubstring(meshFilename, "<domain>",  std::string(domainString));
      ReplaceSubstring(meshFilename, "<time>",    std::string(timeString));

      meshWriter.Write(meshFilename.c_str(), 
        distributedElasticMeshes[domainNumber]->volumeMesh.nodes, 
        distributedElasticMeshes[domainNumber]->volumeMesh.cells, std::vector< typename AdditionalCellInfo<Space>:: template AuxInfo<char> >(),
        &(distributedElasticMeshes[domainNumber]->volumeMesh.isCellAvailable));
    }

    if(settings.snapshots[snapshotIndex].contacts.used)
    {
      std::string contactsFilename = settings.snapshots[snapshotIndex].contacts.filename;
      ReplaceSubstring(contactsFilename, "<step>",    std::string(stepString));
      ReplaceSubstring(contactsFilename, "<domain>",  std::string(domainString));
      ReplaceSubstring(contactsFilename, "<time>",    std::string(timeString));
      contactsWriter.Write(contactsFilename,
        distributedElasticMeshes[domainNumber], 
        settings.snapshots[snapshotIndex].contacts.faceTypesToDraw, 
        ContactVtkWriter<Space, FunctionSpace>::Contact);
    }

    if (settings.snapshots[snapshotIndex].boundaries.used)
    {
      std::string boundariesFilename = settings.snapshots[snapshotIndex].boundaries.fileName;
      ReplaceSubstring(boundariesFilename, "<step>", std::string(stepString));
      ReplaceSubstring(boundariesFilename, "<domain>", std::string(domainString));
      ReplaceSubstring(boundariesFilename, "<time>", std::string(timeString));

      contactsWriter.Write(boundariesFilename,
        distributedElasticMeshes[domainNumber], 
        settings.snapshots[snapshotIndex].boundaries.faceTypesToDraw,
        ContactVtkWriter<Space, FunctionSpace>::Boundary);
    }

    if (settings.snapshots[snapshotIndex].cellInfos.used)
    {
      std::string cellInfosFilename = settings.snapshots[snapshotIndex].cellInfos.fileName;
      ReplaceSubstring(cellInfosFilename, "<step>", std::string(stepString));
      ReplaceSubstring(cellInfosFilename, "<domain>", std::string(domainString));
      ReplaceSubstring(cellInfosFilename, "<time>", std::string(timeString));

      cellInfosWriter.Write(cellInfosFilename, distributedElasticMeshes[domainNumber]);
    }
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::AnalyzeSnapshotData(Scalar currTime, IndexType snapshotIndex, IndexType stepIndex, Overload<Space3>)
{
  printf("This is debug processing code. Don't forget to turn it off for regular tasks.\n");
}

template<typename Space, unsigned int order>
void Task<Space, order>::AnalyzeSnapshotData(Scalar currTime, IndexType snapshotIndex, IndexType stepIndex, Overload<Space2>)
{
  printf("This is debug processing code. Don't forget to turn it off for regular tasks.\n");

  IniStateMaker<ElasticSpace> *dstStateMaker = stateMakers[0]; //state to match to

  int sampleWidth = 512; //integration area resolution
  int sampleHeight = 512;

  Scalar eps = Scalar(1e-3); //to compute integral only in the mesh and avoid interpolating outside of it
  Vector boxPoint1 = Vector(eps, eps); //integration area
  Vector boxPoint2 = Vector(Scalar(1) - eps, Scalar(1) - eps);

  std::string analysisFileName = "out/analysis<step>.txt"; //analysis file name

  if(GetCurrentNodeDomainsCount() != 1) return;

  char stepString[256];
  char timeString[256];

  // int domainIndex = 0;
  sprintf(stepString, "%.6d", static_cast<int>(stepIndex));
  sprintf(timeString, "%.5f", currTime);

  ReplaceSubstring(analysisFileName, "<step>",    std::string(stepString));
  ReplaceSubstring(analysisFileName, "<time>",    std::string(timeString));

  std::fstream analysisDataFile;
  analysisDataFile.open(analysisFileName.c_str(), std::fstream::out);


  std::vector<Elastic> sampleData;
  sampleData.resize(sampleWidth * sampleHeight);


  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    // snapshot
    // AABB2 snapshotArea = distributedElasticMeshes[domainNumber]->volumeMesh.GetComputationArea();
    distributedElasticMeshes[domainNumber]->MakeSnapshot(sampleData.data(), 
      sampleWidth, 
      sampleHeight, 
      boxPoint1, 
      boxPoint2, false);

    //distributedElasticMeshes[domainNumber]->
    Vector2 stepSize = Vector((boxPoint2.x - boxPoint1.x) / Scalar(sampleWidth - 1), (boxPoint2.y - boxPoint1.y) / Scalar(sampleHeight - 1));

    Scalar l2DiffSqr = 0;
    Scalar lInfDiff = 0;
    for(int y = 0; y < sampleHeight; y++)
    {
      {
        for(int x = 0; x < sampleWidth; x++)
        {
          Vector point(boxPoint1.x + Scalar(x) * stepSize.x,
                       boxPoint1.y + Scalar(y) * stepSize.y);
          assert(stateMakers.size() > 0);

          typename ElasticSpace::MediumParametersType material;
          material.SetFromVelocities(2, 1, 1);

          assert(stateMakers.size() > 0);
          Elastic dstValue = dstStateMaker->GetValue(point, material.lambda, material.mju, material.invRho);
          Elastic value = sampleData[x + y * sampleWidth];

          Scalar err = (dstValue.GetVelocity() - value.GetVelocity()).Len();
          l2DiffSqr += (stepSize.x * stepSize.y) * Sqr(err);
          lInfDiff = std::max(lInfDiff, abs(err));
        }
      }
    }
    Scalar l2Diff = sqrt(l2DiffSqr);
    std::stringstream analysisLineStr;
    analysisLineStr << 
      "L2 error: " << l2Diff << std::endl << 
      "L Inf diff: " << lInfDiff << std::endl <<
      "Domain: " << domainNumber << std::endl << 
      "Step: " << stepIndex << std::endl << 
      "Time " << currTime << std::endl;

    analysisDataFile << analysisLineStr.str();
  }
  analysisDataFile.close();
}

template<typename Space, unsigned int order>
void Task<Space, order>::LoadGeom(IndexType domainNumber)
{
  LoadGeom(Overload<Space>(), domainNumber);
}

template<typename Space, unsigned int order>
void Task<Space, order>::LoadMeshes()
{
  assert(Space::Dimension == settings.configDimsCount);

  meshes.resize(GetCurrentNodeDomainsCount());
  syncData.resize(GetCurrentNodeDomainsCount());
  isContactBroken.resize(GetCurrentNodeDomainsCount());

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {

    meshes[domainNumber] = new DistributedMeshIO<Space>(GetDomainsCount());
    char domainString[256];
    sprintf(domainString, "%.3d", static_cast<int>(nodesSchedule[GetNodeId()].domainsIndices[domainNumber]));

    std::string meshName = AddExtensionToFileName(settings.mesh.meshFileName, ".mesh");
    ReplaceSubstring(meshName, "<domain>", domainString);

    printf("Loading %s geom: node %d; domain %d from file %s\n",
      (settings.configDimsCount == 2) ? "2d" : "3d",
      static_cast<int>(GetNodeId()),
      static_cast<int>(nodesSchedule[GetNodeId()].domainsIndices[domainNumber]),
      meshName.c_str());

    meshes[domainNumber]->Load(meshName);



    IndexType meshCellsCount = meshes[domainNumber]->GetCellsCount();
    cellMediumParams.resize(meshCellsCount);
    internalContactTypes.resize(meshCellsCount);
    typedef typename MeshSettings<Space>::MediumParamsSection MediumParamsSection;
    typedef typename MediumParamsSection::ParamsDescription ParamsDescription;
    ParamsDescription paramsDesc = settings.mesh.mediumParamsSection.paramsDescription;
    switch (paramsDesc.type)
    {
      case ParamsDescription::PerSubmesh:
      {
        typename MediumParamsSection::PerSubmeshInfo perSubmeshInfo =
          settings.mesh.mediumParamsSection.perSubmeshInfos[paramsDesc.infoIndex];

        std::string paramsFileName = perSubmeshInfo.fileName;
        if(paramsFileName == "")
        {
          paramsFileName = settings.mesh.meshFileName;
        }
        paramsFileName = AddExtensionToFileName(paramsFileName, ".params");

        ReplaceSubstring(paramsFileName, "<domain>", domainString);

        FILE* paramsFile = fopen(paramsFileName.c_str(), "rb");
        if (paramsFile)
        {
          for(IndexType cellIndex = 0; cellIndex < meshCellsCount; ++cellIndex)
          {
            char currCellSubmeshIndex = char(-1);
            IndexType bytesRead = fread(&currCellSubmeshIndex, sizeof(char), 1, paramsFile);
            assert(bytesRead == 1);

            for(IndexType submeshNumber = 0; submeshNumber < perSubmeshInfo.submeshParams.size(); ++submeshNumber)
            {
              if(perSubmeshInfo.submeshParams[submeshNumber].submeshIndex == currCellSubmeshIndex)
              {
                cellMediumParams[cellIndex] = MakeElasticMediumParams(perSubmeshInfo.submeshParams[submeshNumber].params);
                internalContactTypes[cellIndex] = perSubmeshInfo.submeshParams[submeshNumber].internalContactType;
              }
            }
          }
          fclose(paramsFile);
        } else
        {
          printf("Can`t open %s file: %s\n", paramsFileName.c_str(), strerror(errno));
          assert(0);
        }
      } break;

      case ParamsDescription::Uniform:
      {
        typename MediumParamsSection::UniformInfo uniformInfo =
          settings.mesh.mediumParamsSection.uniformInfos[paramsDesc.infoIndex];

        for(IndexType cellIndex = 0; cellIndex < meshCellsCount; ++cellIndex)
        {
          cellMediumParams[cellIndex] = MakeElasticMediumParams(uniformInfo.params);
          internalContactTypes[cellIndex] = uniformInfo.internalContactType;
        }
      } break;

      case ParamsDescription::PerCell:
      {
        typename MediumParamsSection::PerCellInfo perCellInfo =
          settings.mesh.mediumParamsSection.perCellInfos[paramsDesc.infoIndex];

        std::string paramsFileName = AddExtensionToFileName(perCellInfo.fileName, ".params");
        
        if(paramsFileName == "")
        {
          paramsFileName = AddExtensionToFileName(settings.mesh.meshFileName, ".params");
        }
        ReplaceSubstring(paramsFileName, "<domain>", domainString);

        FILE *paramsFile = fopen(paramsFileName.c_str(), "rb");

        for(IndexType cellIndex = 0; cellIndex < meshCellsCount; ++cellIndex)
        {
          IndexType currCellSubmeshIndex = IndexType(-1);

          typename MediumParamsSection::MediumParams params;
          IndexType bytesRead = fread(&params, sizeof(params), 1, paramsFile);
          assert(bytesRead == 1);

          cellMediumParams[cellIndex] = MakeElasticMediumParams(params);
          internalContactTypes[cellIndex] = IndexType(0); //TODO: left only for binary compatibility with currently existing per cell info files. Should be changed to from file as well.
        }

        fclose(paramsFile);
      } break;
    }

    for (IndexType vertexIndex = 0; vertexIndex < meshes[domainNumber]->vertices.size(); ++vertexIndex)
    {
      meshes[domainNumber]->vertices[vertexIndex] /= settings.solver.velocityDimensionlessMult;
    }

    LoadGeom(domainNumber);

    for (IndexType cellIndex = 0; cellIndex < distributedElasticMeshes[domainNumber]->volumeMesh.cells.size(); ++cellIndex)
    {
      distributedElasticMeshes[domainNumber]->volumeMesh.cellMediumParameters[cellIndex].initialVolume =
        distributedElasticMeshes[domainNumber]->volumeMesh.GetVolume(cellIndex);
    }

    distributedElasticMeshes[domainNumber]->Initialize(
      settings.solver.allowMovement,
      settings.mesh.collisionWidth,
      settings.mesh.unfoldIterationsCount,
      settings.mesh.minGridHeight,
      settings.solver.tensionErrorMult,
      settings.solver.velocityErrorMult,
      settings.solver.positionErrorMult,
      settings.solver.tensionDimensionlessMult,
      settings.solver.velocityDimensionlessMult);

    SetTransitionInfo(domainNumber, meshes[domainNumber]);

    for (IndexType detectorIndex = 0; detectorIndex < meshes[domainNumber]->detectorsPositions.size(); ++detectorIndex)
    {
      meshes[domainNumber]->detectorsPositions[detectorIndex] /= settings.solver.velocityDimensionlessMult;
    }

    distributedElasticMeshes[domainNumber]->SetDetectors(meshes[domainNumber]->detectorsPositions);
    detectorsData[domainNumber].resize(meshes[domainNumber]->detectorsPositions.size());
    
    isContactBroken[domainNumber].resize(meshes[domainNumber]->GetContactsCount(), false);
  }

  // initialize sync data
  ComputePacketsToReceiveCount();

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    syncData[domainNumber].resize(GetDomainsCount());
    for (IndexType dstDomainIndex = 0; dstDomainIndex < GetDomainsCount(); ++dstDomainIndex)
    {
      syncData[domainNumber][dstDomainIndex].sendingData =
        new char[distributedElasticMeshes[domainNumber]->GetSyncDataSize(dstDomainIndex)];
    }
  }

  BuildBoundaryDescriptions();
  BuildContactDescriptions();

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    if (settings.solver.allowMovement)
    {
      distributedElasticMeshes[domainNumber]->volumeMesh.BuildAABBTree(settings.solver.dynamicContactBox.boxPoint1, settings.solver.dynamicContactBox.boxPoint2);
    }
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::AllocateSnapshotData()
{
  IndexType snapshotSize = FindSnapshotSize();
  snapshotData = new Elastic[snapshotSize];
}

template<typename Space, unsigned int order>
void Task<Space, order>::BuildIniStateMakers()
{
  for (IndexType stateIndex = 0; stateIndex < settings.task.iniStates.size(); ++stateIndex)
  {
    IndexType infoIndex = settings.task.iniStates[stateIndex].infoIndex;
    switch (settings.task.iniStates[stateIndex].type)
    {
      case TaskSettings<Space>::IniState::Box:
      {
        typename TaskSettings<Space>::BoxStateInfo boxInfo = settings.task.boxStateInfos[infoIndex];
        stateMakers.push_back(new BoxIniStateMaker<ElasticSpace>(boxInfo.velocity, AABB(boxInfo.boxPoint1, boxInfo.boxPoint2)));
      }break;
      case TaskSettings<Space>::IniState::Berlage:
      {
        typename TaskSettings<Space>::BerlageStateInfo berlageInfo = settings.task.berlageStateInfos[infoIndex];
        stateMakers.push_back(new BerlageIniStateMaker<ElasticSpace>(berlageInfo.point, berlageInfo.waveVector,
          berlageInfo.waveLength, Scalar(1000), berlageInfo.shear));
      }break;
      case TaskSettings<Space>::IniState::AcousticPlaneWave:
      {
        typename TaskSettings<Space>::AcousticPlaneWaveStateInfo acousticPlaneWaveInfo =
          settings.task.acousticPlaneWaveStateInfos[infoIndex];
        stateMakers.push_back(new AcousticPlaneWaveIniStateMaker<ElasticSpace>(
          acousticPlaneWaveInfo.point, acousticPlaneWaveInfo.waveVector, acousticPlaneWaveInfo.waveLength));
      }break;
      case TaskSettings<Space>::IniState::Ricker:
      {
        typename TaskSettings<Space>::RickerStateInfo rickerInfo = settings.task.rickerStateInfos[infoIndex];
        stateMakers.push_back(new RickerIniStateMaker<ElasticSpace>(rickerInfo.point, rickerInfo.waveVector, rickerInfo.waveLength));
      }break;
      case TaskSettings<Space>::IniState::Radial:
      {
        typename TaskSettings<Space>::RadialStateInfo radialInfo = settings.task.radialStateInfos[infoIndex];
        stateMakers.push_back(new RadialIniStateMaker<ElasticSpace>(radialInfo.point, radialInfo.radialVelocity,
          radialInfo.waveLength, Scalar(1000)));
      }break;
      case TaskSettings<Space>::IniState::Ball:
      {
        typename TaskSettings<Space>::BallStateInfo ballInfo = settings.task.ballStateInfos[infoIndex];
        stateMakers.push_back(new BallIniStateMaker<ElasticSpace>(ballInfo.point, ballInfo.radialVelocity,
          ballInfo.pressure, ballInfo.waveLength, Scalar(1000)));
      }break;
      case TaskSettings<Space>::IniState::Bagel:
      {
        typename TaskSettings<Space>::BagelStateInfo bagelInfo = settings.task.bagelStateInfos[infoIndex];
        stateMakers.push_back(new BagelIniStateMaker<ElasticSpace>(bagelInfo.r1, bagelInfo.r2, bagelInfo.point, bagelInfo.magnitude));
      }break;
      case TaskSettings<Space>::IniState::DirectionalExplosion:
      {
        typename TaskSettings<Space>::DirectionalExplosionStateInfo explosionInfo =
          settings.task.directionalExplosionStateInfos[infoIndex];
        stateMakers.push_back(new DirectionalExplosionIniStateMaker<ElasticSpace>(explosionInfo.point, explosionInfo.velocity,
          explosionInfo.waveLength, Scalar(1000)));
      }break;

      case TaskSettings<Space>::IniState::ConstantPlaneWave:
      {
        typename TaskSettings<Space>::ConstantPlaneWaveStateInfo waveInfo =
          settings.task.constantPlaneWaveStateInfos[infoIndex];
        stateMakers.push_back(new ConstantPlaneWaveIniStateMaker<ElasticSpace>(waveInfo.velocity, waveInfo.center,
          waveInfo.waveLength, waveInfo.shear));
      }break;
    }
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::BuildSourceTerms()
{
  SourceFunctor< Space >* sourceFunctor = nullptr;

  //actually only the first source is used. still looks cool.
  for(IndexType sourceIndex = 0; sourceIndex < settings.task.sourceTerms.size(); ++sourceIndex)
  {
    IndexType infoIndex = settings.task.sourceTerms[sourceIndex].infoIndex;
    switch(settings.task.sourceTerms[sourceIndex].type)
    {
      case TaskSettings<Space>::SourceTerm::ConstantAcceleration:
      {
        typename TaskSettings<Space>::ConstantAccelerationSourceInfo info = settings.task.constantAccelerationSourceInfos[infoIndex];
        sourceFunctor = new ConstantAccelerationFunctor< ElasticSystem<Space> >(info.acceleration, 
          settings.solver.tensionDimensionlessMult, settings.solver.velocityDimensionlessMult);
      }break;
    }
    break;
  }

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    ElasticSystem<Space> *system = distributedElasticMeshes[domainNumber]->GetSystem();
    system->SetSourceFunctor(sourceFunctor);
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::BuildPointSources()
{
  std::vector< PointSource<Space>* > pointSources;

  for (IndexType sourceIndex = 0; sourceIndex < settings.task.pointSources.size(); ++sourceIndex)
  {
    IndexType infoIndex = settings.task.pointSources[sourceIndex].infoIndex;
    switch (settings.task.pointSources[sourceIndex].type)
    {
      case TaskSettings<Space>::PointSource::Force:
      {
        typename TaskSettings<Space>::ForcePointSourceInfo info = settings.task.forcePointSourceInfos[infoIndex];
        PointSource<Space>* source = new ForcePointSource< ElasticSystem<Space> >(info.point,
          info.peakFrequency, info.acceleration, info.latency,
          settings.solver.tensionDimensionlessMult, settings.solver.velocityDimensionlessMult);
        pointSources.push_back(source);
      } break;
      case TaskSettings<Space>::PointSource::Monopole:
      {
        typename TaskSettings<Space>::MonopoleSourceInfo info = settings.task.monopoleSourceInfos[infoIndex];
        PointSource<Space>* source = new MonopoleSource< ElasticSystem<Space> >(info.point,
          info.pressure, info.peakFrequency, info.latency,
          settings.solver.tensionDimensionlessMult, settings.solver.velocityDimensionlessMult);
        pointSources.push_back(source);
      } break;
    } break;
  }

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    ElasticSystem<Space> *system = distributedElasticMeshes[domainNumber]->GetSystem();
    system->SetPointSources(pointSources);
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::LoadInitialState()
{
  if (stateMakers.size() == 0)
  {
    printf("State maker was not created\n");
    return;
  }

  printf("Loading initial state\n");

  for (IndexType stateMakerIndex = 0; stateMakerIndex < stateMakers.size(); ++stateMakerIndex)
  {
    stateMakers[stateMakerIndex]->MakeParamsDimensionless(
      settings.solver.tensionDimensionlessMult, 
      settings.solver.velocityDimensionlessMult);
  }

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    Scalar mult(1.0);
    distributedElasticMeshes[domainNumber]->LoadState(stateMakers, mult);
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::PrintSolverState(const SolverState& solverState, Scalar currTime, Scalar timeStep)
{
  if (solverState.hierarchyLevelsCount > 1)
  {
    printf("step %d, level %d, phase %d, currTime %f, dt = %.8f\n", 
      solverState.globalStepIndex, solverState.hierarchyLevel, solverState.hierarchyPhase, currTime, timeStep);
  } else
  {
    printf("step %d, currTime %f, dt = %.8f\n", 
      solverState.globalStepIndex, currTime, timeStep);
  }
}

template<typename Space, unsigned int order>
Task<Space, order>::Task(): NotifyListener(), ReceiveListener(),
  distributedElasticMeshes(0), network(new NetworkInterface()), packetsToReceive(0), snapshotData(nullptr)
{
}

template<typename Space, unsigned int order>
Task<Space, order>::~Task()
{
  delete network;
  delete [] snapshotData;
  for (IndexType meshIndex = 0; meshIndex < distributedElasticMeshes.size(); ++meshIndex)
  {
    delete distributedElasticMeshes[meshIndex];
    delete meshes[meshIndex];
  }

  for (IndexType srcDomainIndex = 0; srcDomainIndex < syncData.size(); ++srcDomainIndex)
  {
    for (IndexType dstDomainIndex = 0; dstDomainIndex < syncData[srcDomainIndex].size(); ++dstDomainIndex)
    {
      delete [] syncData[srcDomainIndex][dstDomainIndex].sendingData;
    }
  }

  for (IndexType stateMakerIndex = 0; stateMakerIndex < stateMakers.size(); ++stateMakerIndex)
  {
    delete stateMakers[stateMakerIndex];
  }
  profilingDataFile.close();
}

template<typename Space, unsigned int order>
void Task<Space, order>::SaveDetectorsDataToFile()
{
  if (!settings.detectors.used) return;

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    IndexType domainIndex = nodesSchedule[GetNodeId()].domainsIndices[domainNumber];

    char domainString[256];
    sprintf(domainString, "%.3d", static_cast<int>(domainIndex));

    for (IndexType detectorIndex = 0; detectorIndex < distributedElasticMeshes[domainNumber]->detectorsPositions.size(); ++detectorIndex)
    {
      std::string detectorFilename = settings.detectors.filename;
      ReplaceSubstring(detectorFilename, "<domain>", std::string(domainString));

      char detectorString[256];
      sprintf(detectorString, "%.3d", static_cast<int>(detectorIndex));
      ReplaceSubstring(detectorFilename, "<detector>", std::string(detectorString));

      ReplaceSubstring(detectorFilename, "<coords>",
        // transform to dimension coords
        (meshes[domainNumber]->detectorsPositions[detectorIndex] *
        settings.solver.velocityDimensionlessMult).ToString());

      std::fstream seismogramFile;
      seismogramFile.open(detectorFilename.c_str(), std::ios_base::out | std::ios_base::app);
      assert(seismogramFile.is_open());

      for (IndexType sampleIndex = 0; sampleIndex < detectorsData[domainNumber][detectorIndex].size(); ++sampleIndex)
      {
        const Sample& sample = detectorsData[domainNumber][detectorIndex][sampleIndex];
        seismogramFile << sample.time << ";";
        for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          seismogramFile << sample.elastic.values[valueIndex] << ";";
        }
        seismogramFile << std::endl;
      }
      detectorsData[domainNumber][detectorIndex].clear();
      seismogramFile.close();
    }
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::SaveDetectorsData(Scalar currTime, IndexType stepIndex)
{
  if (!settings.detectors.used || stepIndex % settings.detectors.samplesPeriod != 0) return;

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    for (IndexType detectorIndex = 0; detectorIndex < distributedElasticMeshes[domainNumber]->detectorsPositions.size(); ++detectorIndex)
    {
      Sample sample;
      sample.time = currTime;   
      distributedElasticMeshes[domainNumber]->GetDetectorsData(detectorIndex, sample.elastic.values);
      detectorsData[domainNumber][detectorIndex].push_back(sample);
    }
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::UpdateMeshData(char* data)
{
  /*
    packet: dstDomainIndex, data
  */
  IndexType dstDomainIndex = ((IndexType*)data)[0];
  if (solverState.IsUseful())
  {
    distributedElasticMeshes[domainsNumbers[dstDomainIndex]]->UpdateDomainData(data + sizeof(IndexType));
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::OnNotify()
{
  --notificationsToReceive;
}

template<typename Space, unsigned int order>
void Task<Space, order>::OnDataReceive(void *data, int, int)
{
  UpdateMeshData(static_cast<char*>(data));
  ++receivedPackets; 
}

template<typename Space, unsigned int order>
void Task<Space, order>::SynchronizeMeshes()
{
  network->initializeSync();
  double beginSync = MPI_Wtime();

  notificationsToReceive = receivedPackets = 0;

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    IndexType domainIndex = nodesSchedule[GetNodeId()].domainsIndices[domainNumber];
    for (IndexType dstDomainIndex = 0; dstDomainIndex < GetDomainsCount(); ++dstDomainIndex)
    if (!distributedElasticMeshes[domainNumber]->IsSyncDataEmpty(dstDomainIndex))
    {
      distributedElasticMeshes[domainNumber]->BuildSyncData(dstDomainIndex, syncData[domainNumber][dstDomainIndex].sendingData);

      IndexType dstHost = domainsLocation[dstDomainIndex];
      network->sendDataAsync(syncData[domainNumber][dstDomainIndex].sendingData, 
                            static_cast<int>(distributedElasticMeshes[domainNumber]->GetSyncDataSize(dstDomainIndex)), 
                            static_cast<int>(dstHost), this);
       
      ++notificationsToReceive;
    }
  } 
  double endSending = MPI_Wtime();
  SaveProfilingData("   Sending data", beginSync, endSending);
    
  while(receivedPackets != packetsToReceive || notificationsToReceive > 0)
  {
    network->update();
  }
  network->barrier();
  network->waitAll();

  double endSync = MPI_Wtime();
  SaveProfilingData(" Synchronization", beginSync, endSync);
}

template<typename Space, unsigned int order>
void Task<Space, order>::LoadNodesSchedule()
{
  nodesSchedule = settings.schedule.nodesSchedule;
  domainsLocation.resize(GetDomainsCount());

  for (IndexType nodeIndex = 0; nodeIndex < nodesSchedule.size(); ++nodeIndex)
  {
    for (IndexType domainNumber = 0; domainNumber < nodesSchedule[nodeIndex].domainsIndices.size(); ++domainNumber)
    {
      IndexType domainIndex = nodesSchedule[nodeIndex].domainsIndices[domainNumber];
      domainsLocation[domainIndex] = nodeIndex;
    }
  }

  domainsNumbers.resize(GetDomainsCount(), IndexType(-1));
  for (IndexType domainNumber = 0; domainNumber < nodesSchedule[GetNodeId()].domainsIndices.size(); ++domainNumber)
  {
    IndexType domainIndex = nodesSchedule[GetNodeId()].domainsIndices[domainNumber];
    domainsNumbers[domainIndex] = domainNumber;
  }
}

template<typename Space, unsigned int order>
typename Task<Space, order>::IndexType Task<Space, order>::GetDomainsCount() const
{
  return settings.schedule.domainsCount;
}

template<typename Space, unsigned int order>
typename Task<Space, order>::IndexType Task<Space, order>::GetNodeId() const
{
  return network->getID();
}

template<typename Space, unsigned int order>
typename Task<Space, order>::IndexType Task<Space, order>::GetCurrentNodeDomainsCount() const
{
  return nodesSchedule[GetNodeId()].domainsIndices.size();
}

template<typename Space, unsigned int order>
void Task<Space, order>::SetThreadsCount()
{
  #ifdef SINGLE_THREAD
    omp_set_num_threads(1);
    return;
  #endif

  int size, rank;
  int nthreads;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Environment variable OMP_NUM_THREADS is defined only on the 'Mother Superior' node.
  // So we need to send it value to other nodes and apply there.
  if (rank == 0)
  {
    #pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }
    for (int i = 1; i < size; i++) 
    {
      MPI_Send(&nthreads, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  } else 
  {
    MPI_Recv(&nthreads, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    omp_set_num_threads(nthreads);
  }
  #pragma omp parallel private(nthreads)
  {
    nthreads = omp_get_num_threads();
  }
}

template<typename Space, unsigned int order>
void Task<Space, order>::ComputePacketsToReceiveCount()
{
  struct CacheNotifyListener: public NotifyListener
  {
    virtual ~CacheNotifyListener() = default;
    virtual void OnNotify() override
    {
      --notificationsRemained;
    }
    IndexType notificationsRemained;
  } cacheNotifyListener;

  struct CacheReceiveListener: public ReceiveListener
  {
    virtual ~CacheReceiveListener() = default;
    virtual void OnDataReceive(void *data, int size, int hostId) override
    {
      ++packetsToReceive;
    }
    IndexType packetsToReceive;
  } cacheReceiveListener;

  network->setReceiveListener(&cacheReceiveListener);
  cacheNotifyListener.notificationsRemained = 0;
  cacheReceiveListener.packetsToReceive = 0;
  network->initializeSync();

  char dummyData;
  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    for (IndexType dstDomainIndex = 0; dstDomainIndex < GetDomainsCount(); ++dstDomainIndex)
    {
      char nonEmpty = !distributedElasticMeshes[domainNumber]->IsSyncDataEmpty(dstDomainIndex);
      IndexType dstHost = domainsLocation[dstDomainIndex];
      if (nonEmpty)
      {
        network->sendDataAsync(&dummyData, sizeof(char), static_cast<int>(dstHost), &cacheNotifyListener);
        ++cacheNotifyListener.notificationsRemained;
      }
    }
  }

  bool notified = false;
  while (!network->checkAllReadyAsync())
  {
    network->update();
    if (cacheNotifyListener.notificationsRemained == 0 && !notified)
    {
      network->notifyReady();
      notified = true;
    }
  }
  network->barrier();
  network->waitAll();

  network->setReceiveListener(this);
  packetsToReceive = cacheReceiveListener.packetsToReceive;
}

template<typename Space, unsigned int order>
void Task<Space, order>::SaveProfilingData(const std::string& info, double begin, double end)
{
  #ifdef PROFILING
    double elapsed_msecs = end - begin;
    profilingDataFile << '[' << GetNodeId() << "," << solverState.globalStepIndex << "] " << info << ": " << elapsed_msecs << std::endl;
  #endif
}

template<typename Space, unsigned int order>
void Task<Space, order>::WriteEnergyAndImpulseDeviation(const std::vector<Scalar>& initialEnergies,
  const std::vector<Vector>& initialImpulses, Scalar currTime)
{
  #ifdef WRITE_ENERGY_AND_IMPULSE
  std::fstream file("out/energy.txt", std::fstream::out | std::fstream::app);
  if (!file.is_open()) return;
  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    Vector totalImpulse = distributedElasticMeshes[domainNumber]->GetTotalImpulse();

    Scalar kineticEnergy, potentialEnergy;
    Scalar totalEnergy = distributedElasticMeshes[domainNumber]->GetTotalEnergy(kineticEnergy, potentialEnergy);
    file << currTime << " " 
                     << "Energy deviation: " << totalEnergy - initialEnergies[domainNumber] << "; "
                     << "Energy: " << totalEnergy << "; " 
                     << "Impulse deviation: " << (totalImpulse - initialImpulses[domainNumber]).Len() << "; "
                     << "Impulse: ";
    for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
    {
      file << totalImpulse[dimIndex] << " ";
    }
    
    file << "Mass: " << distributedElasticMeshes[domainNumber]->GetTotalMass();
    file << " Kinetic E: " << kineticEnergy << ";";
    file << " Potential E: " << potentialEnergy;
    file << std::endl;
  }
  file.close();
  #endif // WRITE_ENERGY_AND_IMPULSE
}


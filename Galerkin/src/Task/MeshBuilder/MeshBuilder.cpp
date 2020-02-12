#include "../../Maths/Spaces.h"
#include "../GeomMesh/GeomMesh/GeomMesh.h"
#include "SpecificMeshBuilders/SpecificMeshBuilders.h"
#include "../GeomMesh/TimeHierarchyLevelsManager.h"

#include "../GeomMesh/MeshIO/Distributed/DistributedMeshIO.h"
#include "../GeomMesh/MeshIO/Local/MeshIO.h"
#include "MeshSplitter/MeshSplitter.h"
#include "MeshDistributor.h"
#include "MeshChecker.h"

#include <sstream>
#include <string>
#include <vector>
#include "../../IO/Vtk/MeshVtkWriter.h"
#include "../Task/SettingsParser/SettingsParser.h"
#include "../../DifferentialSolvers/SolversFactory.h"

/* don't change. meshbuilder compilation is fast enough and it's extremely 
difficult to detect Space2/Space3 mismatch because nothing crashes */

#define SPACE_FROM_SETTINGS 
//typedef Space2 DefaultSpace;

template<typename Space>
struct MeshBuilderTask
{
public:
  void Run()
  {
    settings.Parse("task.xml");

    std::string meshFileName = settings.mesh.meshFileName;
    if(settings.meshBuilder.meshFileName != "")
    {
      meshFileName = settings.meshBuilder.meshFileName;
    }
    this->meshBaseName = meshFileName;

    assert(settings.configDimsCount == Space::Dimension);


    // build geom
    mesh = new MeshIO<Space>();
    geomMesh = new GeomMesh<Space>();
    nonDuplicatedVerticesMesh = new MeshIO<Space>();
    BuildGeom();

    meshSplitter    = new MeshSplitter<Space>(mesh, geomMesh);

    BuildGlobalMediumParams(geomMesh);

    // distribute
    distributor = new MetisDistributor<Space>(); 
    // destributor = new SimpleRectDistributor<Space2>(); 
    Distribute(domains);

    IndexType paramsPerCellSize = globalMediumParams.size() / mesh->GetCellsCount();
    SaveDomains(paramsPerCellSize);
    delete meshSplitter;
    delete mesh;
    delete geomMesh;
    delete nonDuplicatedVerticesMesh;

    // mesh checking
    MeshChecker<Space> meshChecker(meshBaseName, settings.schedule.domainsCount);
    meshChecker.Check();
    printf("\nDone");
  }
private:
  typedef typename Space::IndexType                              IndexType;
  typedef typename Space::Scalar                                 Scalar;
  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<Scalar>    ScalarCellInfo;
  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<IndexType> IndexCellInfo;
  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<char>      CharCellInfo;

  Settings<Space>                         settings;
  MeshIO<Space>*                          mesh;
  GeomMesh<Space>*                        geomMesh;
  MeshIO<Space>*                          nonDuplicatedVerticesMesh;
  BasicMeshBuilder<Space>*                builder;
  MeshDestributor<Space>*                 distributor;
  MeshSplitter<Space>*                    meshSplitter;
  std::vector< DistributedMeshIO<Space> > domains;
  std::vector<char>                       globalMediumParams;

  std::string meshBaseName;

  void BuildGeom()
  {
    std::cout << "Building mesh...\n";

    bool duplicateNodes = settings.solver.allowDiscreteDestruction || settings.solver.allowContinuousDestruction;
    builder = new SalomeMeshBuilder<Space>(settings.meshBuilder.salomeFileName, duplicateNodes, Scalar(0.0));
    builder->BuildMesh(mesh, geomMesh, nonDuplicatedVerticesMesh);

    std::string globalMeshName = meshBaseName;
    ReplaceSubstring(globalMeshName, "<domain>", "");
    mesh->Save(AddExtensionToFileName(globalMeshName, ".mesh"));
  }

  void BuildGlobalMediumParams(const GeomMesh<Space>* geomMesh)
  {
    printf("Writing global mesh to *.vtk file...\n");
    builder->BuildMediumParams(mesh, &globalMediumParams);
    delete builder;


    std::string globalMeshName = meshBaseName;
    ReplaceSubstring(globalMeshName, "<domain>", "");
    mesh->Save(AddExtensionToFileName(globalMeshName, ".mesh"));

    MeshVtkWriter<Space, CharCellInfo> meshWriter;
    if (globalMediumParams.size() == mesh->GetCellsCount())
    {
      meshWriter.Write(globalMeshName + "_globalMesh.vtk", mesh->vertices, mesh->indices, globalMediumParams);
    } else
    {
      meshWriter.Write(globalMeshName + "_globalMesh.vtk", mesh->vertices, mesh->indices);
    }
  }

  void Distribute(std::vector< DistributedMeshIO<Space> >& domains)
  {
    const IndexType domainsCount = settings.schedule.domainsCount;
    const IndexType hierarchyLevelsCount = settings.solver.hierarchyLevelsCount;
    const IndexType solverPhasesCount = SolversFactory<Scalar, IndexType>::Build(settings.solver.integrator)->GetPhasesCount();

    printf("Distributing...\n");
    IndexType lastComputationalCost = std::numeric_limits<IndexType>::max();

    std::vector< IndexType >  cellWeights(mesh->GetCellsCount(), 1);
    std::vector< Scalar >     cellMaxWaveSpeeds(mesh->GetCellsCount(), Scalar(1.0));
  
    IndexType iterationsCount = 0;
    const IndexType MaxIterationsCount = 20;
    bool efficientStep = true;

    while (efficientStep && iterationsCount < MaxIterationsCount)
    {
      efficientStep = false;
      std::vector<IndexType> cellsDomainIds;

      distributor->Distribute(nonDuplicatedVerticesMesh, domainsCount, cellWeights, settings.meshBuilder.partitionAlgorithmName, &cellsDomainIds);
      std::vector<IndexType> updatedCellsDomainIds(cellsDomainIds.size());

      for (IndexType cellIndex = 0; cellIndex < geomMesh->cells.size(); ++cellIndex)
      {
        updatedCellsDomainIds[geomMesh->updatedCellIndices[cellIndex]] = cellsDomainIds[cellIndex];
      }

      TimeHierarchyLevelsManager<Space> timeHierarchyLevelsManager;
      timeHierarchyLevelsManager.Initialize(mesh->GetCellsCount(), hierarchyLevelsCount, solverPhasesCount, false);

      for (IndexType cellIndex = 0; cellIndex < mesh->GetCellsCount(); ++cellIndex)
      {
        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
        {
          IndexType correspondingCellIndex = 
            geomMesh->GetCorrespondingCellIndex(cellIndex, faceNumber);

          if (correspondingCellIndex != IndexType(-1) && updatedCellsDomainIds[cellIndex] != updatedCellsDomainIds[correspondingCellIndex])
          {
            timeHierarchyLevelsManager.SetLevel(cellIndex, 0);
          }
        }
      }

      Scalar minTimeStep = std::numeric_limits<Scalar>::max();
      for (IndexType cellIndex = 0; cellIndex < geomMesh->cells.size(); ++cellIndex)
      {
        Scalar cellTimeStep = geomMesh->GetMinHeight(cellIndex) / cellMaxWaveSpeeds[cellIndex];
        minTimeStep = std::min(minTimeStep, cellTimeStep);
      }

      timeHierarchyLevelsManager.BuildTimeHierarchyLevels(geomMesh, cellMaxWaveSpeeds, settings.solver.allowMovement, minTimeStep);
      IndexType maxComputationalCost = timeHierarchyLevelsManager.GetComputationalTotalCost(updatedCellsDomainIds, domainsCount);

      if (maxComputationalCost < lastComputationalCost)
      {
        efficientStep = true;
        lastComputationalCost = maxComputationalCost;
        for (IndexType cellIndex = 0; cellIndex < mesh->GetCellsCount(); ++cellIndex)
        {
          cellWeights[cellIndex] = timeHierarchyLevelsManager.GetComputationalCost(geomMesh->updatedCellIndices[cellIndex]);
        }
        timeHierarchyLevelsManager.SaveToVtk(geomMesh, iterationsCount);
      }
      ++iterationsCount;
      std::cout << "  iteration - " << iterationsCount << "; cost - " << maxComputationalCost << ";" << std::endl;
    }

    std::cout << "Splitting..." << std::endl;
    std::vector<IndexType> cellsDomainIds;
    distributor->Distribute(nonDuplicatedVerticesMesh, domainsCount, cellWeights, settings.meshBuilder.partitionAlgorithmName, &cellsDomainIds);
    delete distributor;

    MeshVtkWriter<Space, IndexCellInfo> writer;

    std::string domainsName = meshBaseName;
    ReplaceSubstring(domainsName, "<domain>", "");
    writer.Write(domainsName + "_domains.vtk", mesh->vertices, mesh->indices, cellsDomainIds);

    meshSplitter->Split(cellsDomainIds, domainsCount, &domains, settings.solver.allowMovement);

    // compute partitioning quality
    IndexType maxDomainSize = 0;
    IndexType minDomainSize = geomMesh->cells.size();
    for (IndexType domainIndex = 0; domainIndex < domains.size(); ++domainIndex)
    {
      minDomainSize = std::min(minDomainSize, domains[domainIndex].GetCellsCount());
      maxDomainSize = std::max(maxDomainSize, domains[domainIndex].GetCellsCount());
    }
    std::cout << "Partitioning quality: " << Scalar(maxDomainSize - minDomainSize) / minDomainSize << "\n";

    printf("\n\n");
  }

  void SaveDomains(IndexType paramsPerCellSize)
  {
    printf("Saving domains...\n");
    geomMesh->BuildOriginalCellIndices(mesh->indices.data(), mesh->GetCellsCount());

    for (size_t domainIndex = 0; domainIndex < domains.size(); ++domainIndex)
    {
      std::cout << "  Saving domain " << domainIndex << std::endl;
      char domainString[256];
      sprintf(domainString, "%.3lu", domainIndex);

      std::string domainMeshName = meshBaseName;
      ReplaceSubstring(domainMeshName, "<domain>", std::string(domainString));

      domains[domainIndex].Save(AddExtensionToFileName(domainMeshName, ".mesh"));
      domains[domainIndex].SaveContacts(domainMeshName + "_contacts.vtk");
      domains[domainIndex].SaveBoundaries(domainMeshName + "_boundaries.vtk");
      domains[domainIndex].SaveDetectors(domainMeshName + "_detectors.vtk");
      domains[domainIndex].SaveTransitionInfos(domainMeshName + "_transitionInfos.vtk");

      std::vector<char> mediumParams(domains[domainIndex].GetCellsCount() * paramsPerCellSize, 0);
    
      if (paramsPerCellSize > 0)
      {
        for (IndexType cellIndex = 0; cellIndex < meshSplitter->domainsInfos[domainIndex].cellsIndices.globalIndices.size(); ++cellIndex)
        {
          IndexType updatedGlobalCellIndex = meshSplitter->domainsInfos[domainIndex].cellsIndices.globalIndices.at(cellIndex);
          IndexType originalGlobalCellIndex = geomMesh->originalCellIndices[updatedGlobalCellIndex];
        
          std::copy(globalMediumParams.data() +  originalGlobalCellIndex      * paramsPerCellSize,
                    globalMediumParams.data() + (originalGlobalCellIndex + 1) * paramsPerCellSize,
                    mediumParams.data()       + cellIndex                     * paramsPerCellSize);
        }
      }

      if (mediumParams.size() > 0)
      {
        std::string paramsFileName = AddExtensionToFileName(meshBaseName, ".params");
        ReplaceSubstring(paramsFileName, "<domain>", std::string(domainString));
        FILE *paramsFile = fopen(paramsFileName.c_str(), "wb");
        assert(paramsFile != NULL);
        fwrite(mediumParams.data(), 1, mediumParams.size(), paramsFile);
        fclose(paramsFile);
      }

      std::string vtkMeshFilename = meshBaseName + "_geom.vtk";
      ReplaceSubstring(vtkMeshFilename, "<domain>", std::string(domainString));
      MeshVtkWriter<Space> writer;

      std::vector<CharCellInfo> cellInfos(domains[domainIndex].GetCellsCount());
      for (IndexType cellIndex = 0; cellIndex < domains[domainIndex].GetCellsCount(); ++cellIndex)
      {
        cellInfos[cellIndex].count = 1;
        cellInfos[cellIndex].data = &(mediumParams[cellIndex]);
      }

      writer.Write(vtkMeshFilename, domains[domainIndex].vertices, domains[domainIndex].indices, cellInfos);
    }
    std::cout << "\n";
  
    // save detectors locations (global index -> domain index + detector number)
    std::fstream detectorsLocations(settings.detectors.locationsFileName.c_str(), std::ios_base::out);

    assert(detectorsLocations.fail() == false);

    detectorsLocations << meshSplitter->detectorsLocations.size() << "\n";
    for (IndexType detectorIndex = 0; detectorIndex < meshSplitter->detectorsLocations.size(); ++detectorIndex)
    {
      detectorsLocations << meshSplitter->detectorsLocations[detectorIndex].domainIndex << " "
                         << meshSplitter->detectorsLocations[detectorIndex].localIndex << "\n";
    }
    detectorsLocations.close();
  }
};

int main()
{
  BasicSettings basicSettings; //Space2 settings reader should read Space3 settings file just fine. probably.
  basicSettings.Parse("task.xml");

  #ifdef SPACE_FROM_SETTINGS
  {
    switch(basicSettings.configDimsCount)
    {
      case 2:
      {
        MeshBuilderTask<Space2> meshBuilderTask;
        meshBuilderTask.Run();
      }break;
      case 3:
      {
        MeshBuilderTask<Space3> meshBuilderTask;
        meshBuilderTask.Run();
      }break;
      default: std::cerr << "Wrong dims count in config"; break;
    }
  }
  #else
    MeshBuilderTask<DefaultSpace> meshBuilderTask;
    meshBuilderTask.Run();
  #endif
  
  return 0;
}

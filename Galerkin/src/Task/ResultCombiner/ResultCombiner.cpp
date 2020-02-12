#include "../../Utils/Utils.h"
#include "../../IO/Vtk/BasicVtkWriter.h"
#include "../../IO/Vtk/SnapshotVtkWriter.h"
#include "../../IO/Vtk/VtkReader.h"
#include "../../IO/Segy/SegySeismo.h"

#include "../Task/SettingsParser/SettingsParser.h"
#include "../ElasticSystem/ElasticSystem.h"
#include "SampleIO.h"

#include <fstream>
#include <vector>
#include <stdio.h>

//typedef Space3 DefaultSpace;

#define SPACE_FROM_SETTINGS

template <typename Space>
class ResultCombinerTask
{
public:
  SPACE_TYPEDEFS

  struct Location
  {
    IndexType domainIndex;
    IndexType localIndex;
  };

  typedef ElasticSystem<Space>            ElasticSystemType;
  typedef typename ElasticSystemType::ElasticSpace ElasticSpace;
  typedef typename ElasticSpace::Elastic           Elastic;

  const static int dimsCount = ElasticSystemType::dimsCount;

  void Run()
  {
    settings.Parse("task.xml");

    assert(settings.configDimsCount == Space::Dimension);


    IndexType domainsCount = settings.schedule.domainsCount;

    std::vector<Elastic> elastic, result;
    SnapshotVtkWriter<ElasticSpace> snapshotWriter;
    VtkReader<ElasticSpace> snapshotReader;

    for (IndexType snapshotIndex = 0; snapshotIndex < settings.snapshots.size(); ++snapshotIndex)
    {
      IndexType snapshotSize = ::ComputeSnapshotSize<Space>(settings.snapshots[snapshotIndex].data.origin,
        settings.snapshots[snapshotIndex].data.spacing,
        AABB(settings.snapshots[snapshotIndex].data.boxPoint1,
        settings.snapshots[snapshotIndex].data.boxPoint2));
      if (result.size() < snapshotSize) result.resize(snapshotSize);

      for (int stepIndex = 0; stepIndex < settings.resultCombiner.snapshotsCount; ++stepIndex)
      {
        #pragma omp parallel for
        for (int elasticIndex = 0; elasticIndex < int(snapshotSize); ++elasticIndex)
        {
          result[elasticIndex].SetZeroValues();
        }

        char stepString[256];
        sprintf(stepString, "%.6d", stepIndex);
        bool isOk = true;

        for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
        {
          char domainString[256];
          sprintf(domainString, "%.3lu", domainIndex);

          if (settings.snapshots[snapshotIndex].data.used)
          {
            std::string dataFilename = settings.snapshots[snapshotIndex].data.filename;
            ReplaceSubstring(dataFilename, "<step>", std::string(stepString));
            ReplaceSubstring(dataFilename, "<domain>", std::string(domainString));
            std::string pathName = dataFilename;
            bool res = false;

            if (EndWith(dataFilename, ".vti"))
            {
              res = snapshotReader.Read(pathName, &elastic);
            }
            if (!res)
            {
              isOk = false;
              break;
            }

            #pragma omp parallel for
            for (int elasticIndex = 0; elasticIndex < int(snapshotSize); ++elasticIndex)
            {
              for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
              {
                if (fabs(elastic[elasticIndex].values[valueIndex]) > fabs(result[elasticIndex].values[valueIndex]))
                {
                  result[elasticIndex].values[valueIndex] = elastic[elasticIndex].values[valueIndex];
                }
              }
            }
            if (remove(pathName.c_str()) != 0)
            {
              std::cout << "Can`t remove " << pathName << std::endl;
            }
          }
        }

        if (isOk)
        {
          std::string dataFilename = settings.snapshots[snapshotIndex].data.filename;
          ReplaceSubstring(dataFilename, "<step>", std::string(stepString));
          ReplaceSubstring(dataFilename, "<domain>", "");
          ReplaceSubstring(dataFilename, ".vtk", ".vti");
          std::cout << dataFilename << std::endl;

          Vector boxPoint1 = settings.snapshots[snapshotIndex].data.boxPoint1;
          Vector boxPoint2 = settings.snapshots[snapshotIndex].data.boxPoint2;
          Vector origin = settings.snapshots[snapshotIndex].data.origin;
          Vector spacing = settings.snapshots[snapshotIndex].data.spacing;

          snapshotWriter.Write((dataFilename).c_str(),
            &result[0],
            origin, spacing,
            AABB(boxPoint1, boxPoint2),
            AABB(boxPoint1, boxPoint2),
            settings.snapshots[snapshotIndex].data.writeVelocity,
            settings.snapshots[snapshotIndex].data.writeTension,
            settings.snapshots[snapshotIndex].data.writeWaves);
        }
      }
    }
    CombineDetectors();
  }
private:
  Settings<Space> settings;

  void CombineDetectors()
  {
    SPACE_TYPEDEFS

    std::cout << "Combine detectors\n";

    std::string detectorsFileNamePattern = settings.detectors.filename;
    if(settings.resultCombiner.detectorsFileName != "")
      detectorsFileNamePattern = settings.resultCombiner.detectorsFileName;

    std::string locationsFile = settings.resultCombiner.detectorsLocationsFile;
    if(locationsFile == "")
    {
      locationsFile = settings.detectors.locationsFileName;
    }
    std::fstream detectorsLocationsFile(locationsFile.c_str(), std::ios_base::in);
    assert(detectorsLocationsFile.fail() == false);

    IndexType detectorsCount = 0;
    detectorsLocationsFile >> detectorsCount;

    std::vector<Location> detectorsLocations(detectorsCount, Location());

    for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
    {
      detectorsLocationsFile >>
        detectorsLocations[detectorIndex].domainIndex >>
        detectorsLocations[detectorIndex].localIndex;
    }

    std::vector<FILE*> detectorFiles(detectorsCount, NULL);

    for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
    {
      char detectorString[9];
      sprintf(detectorString, "%.3lu", detectorsLocations[detectorIndex].localIndex);
      char domainString[9];
      sprintf(domainString, "%.3lu", detectorsLocations[detectorIndex].domainIndex);

      std::string fileName = detectorsFileNamePattern;
      ReplaceSubstring(fileName, "<detector>", std::string(detectorString));
      ReplaceSubstring(fileName, "<domain>", std::string(domainString));
      if ((detectorFiles[detectorIndex] = fopen(fileName.c_str(), "r")) == NULL)
      {
        std::cerr << "Can`t open file: " << fileName << "\n";
      }
      assert(detectorFiles[detectorIndex] != NULL);
    }

    std::fstream outputVelocity;
    std::fstream outputPressure;

    if(settings.resultCombiner.velocityCsvName != "")
    {
      printf("Making velocity csv file\n");
      outputVelocity.open(settings.resultCombiner.velocityCsvName.c_str(), std::ios_base::out);
      assert(outputVelocity.fail() == false);
    }
    if(settings.resultCombiner.pressureCsvName != "")
    {
      printf("Making pressure csv file\n");
      outputPressure.open(settings.resultCombiner.pressureCsvName.c_str(), std::ios_base::out);
      assert(outputPressure.fail() == false);
    }

    if(outputVelocity.is_open())
      outputVelocity << "Time;";
    if(outputPressure.is_open())
      outputPressure << "Time;";

    for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
    {
      if(outputVelocity.is_open())
        WriteVelocityHandle<Space>(outputVelocity);
      if(outputPressure.is_open())
        outputPressure << "Vx;";
    }
    if(outputVelocity.is_open())
      outputVelocity << "\n";
    if(outputPressure.is_open())
      outputPressure << "\n";

    bool eof;
    do
    {
      eof = false;
      for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
      {
        Vector v = Vector::zero();
        Scalar time = 0;
        Tensor sigma;

        ReadSample<Space>(detectorFiles[detectorIndex], time, v, sigma);

        if (feof(detectorFiles[detectorIndex]))
        {
          eof = true;
          continue;
        }

        if (detectorIndex == 0)
        {
          if(outputVelocity.is_open())
            outputVelocity << time << ";";
          if(outputPressure.is_open())
            outputPressure << time << ";";
        }
        WriteSample<Space>(outputVelocity, outputPressure, time, v, sigma); //has is_open() inside
      }
      if(outputVelocity.is_open())
        outputVelocity << "\n";
      if(outputPressure.is_open())
        outputPressure << "\n";
    } while (!eof);

    if(outputVelocity.is_open())
      outputVelocity.close();
    if(outputPressure.is_open())
      outputPressure.close();
    // convert *.csv to *.segy
    typedef typename CombinedSeismogramm<float, Space::Dimension>::Elastic SeismoElastic;
    CombinedSeismogramm < float, Space::Dimension > combinedSeismogram = CombinedSeismogramm < float, Space::Dimension >(1);

    bool saveSegY = settings.resultCombiner.velocityCoordSegyName != "";
    bool saveSegYDiff = settings.resultCombiner.velocityDiffCoordSegyName != "";
    if(settings.resultCombiner.velocityCsvName != "" && (saveSegY || saveSegYDiff))
    {
      printf("Combining seismograms into segy\n");
      std::vector<std::string> srcScvFiles;
      srcScvFiles.push_back(settings.resultCombiner.velocityCsvName);

      const std::string dimNames = "xyz";
      std::vector<std::string> dstSegyFiles;
      std::vector<std::string> dstSegyDiffFiles;

      if(saveSegY)
      {
        for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
        {
          std::string fileName = settings.resultCombiner.velocityCoordSegyName;
          ReplaceSubstring(fileName, "<coord>", std::string(1, dimNames[dimIndex]));

          dstSegyFiles.push_back(fileName);
        }
      }
      if(saveSegYDiff)
      {
        for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
        {
          std::string fileName = settings.resultCombiner.velocityDiffCoordSegyName;
          ReplaceSubstring(fileName, "<coord>", std::string(1, dimNames[dimIndex]));

          dstSegyDiffFiles.push_back(fileName);
        }
      }


      combinedSeismogram.Load(CSV, srcScvFiles);

      if(saveSegY)
        combinedSeismogram.Save(SEG_Y, dstSegyFiles);

      if(settings.resultCombiner.velocityRefCsvName != "" && saveSegYDiff)
      {

        std::vector<std::string> refPath;
        refPath.push_back(settings.resultCombiner.velocityRefCsvName);
        CombinedSeismogramm < float, Space::Dimension > refSeismogram = CombinedSeismogramm < float, Space::Dimension >(1);
        refSeismogram.Load(CSV, refPath, &combinedSeismogram);
        printf("Subtracting ref seismogram\n");
        refSeismogram -= combinedSeismogram;

        refSeismogram.Save(SEG_Y, dstSegyDiffFiles);
      }

    }
    printf("Done");
  }
};

int main()
{
  BasicSettings basicSettings; //Space2 settings reader should read Space3 settings file just fine. probably.
  basicSettings.Parse("task.xml");

  #ifdef SPACE_FROM_SETTINGS
  {
    switch (basicSettings.configDimsCount)
    {
      case 2:
      {
        ResultCombinerTask<Space2> resultCombinerTask;
        resultCombinerTask.Run();
      } break;
      case 3:
      {
        ResultCombinerTask<Space3> resultCombinerTask;
        resultCombinerTask.Run();
      } break;
      default: 
        std::cerr << "Wrong dims count in config"; 
      break;
    }
  }
  #else
    ResultCombinerTask<DefaultSpace> resultCombinerTask;
    resultCombinerTask.Run();
  #endif

  return 0;
}

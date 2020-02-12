#pragma once

#include "../../../../3rdparty/tinyxml/tinyxml.h"
#include "../../../../3rdparty/tinyxml/tinystr.h"
#include <limits>
#include <sstream>
#include "ParserUtil.h"
#include "MeshSettingsParser.h"
#include "ScheduleSettingsParser.h"
#include "SnapshotSettingsParser.h"
#include "SolverSettingsParser.h"
#include "TaskSettingsParser.h"
#include "MeshBuilderSettingsParser.h"
#include "ResultCombinerSettingsParser.h"

struct BasicSettings
{
  virtual ~BasicSettings() = default;
  //std::string fileName;
  int configDimsCount; //dimsCount would probably conflict with Space::dimsCount

  virtual void Parse(const char* fileName);
  static std::string ParseSettingsFileName(const char* fileName);

protected:
  void ParseSettingsFile(const char* fileName);
  static TiXmlElement* GetSettingXmlElement(TiXmlDocument& taskFile, const char* fileName)
  {
    bool taskLoadOkay = taskFile.LoadFile(fileName);

    if (!taskLoadOkay)
    {
      std::cerr << "Loading " << fileName << " fails with following error: " <<
        std::string(taskFile.ErrorDesc()) <<
        " in row " << taskFile.ErrorRow() << std::endl;
      throw;
    }

    TiXmlElement* settingsElement = taskFile.FirstChildElement("Settings");
    if (!settingsElement)
    {
      std::cerr << "There is no Settings element";
      throw;
    }
    return settingsElement;
  }

  static std::string GetPath(const std::string& fileName)
  {
    std::string path;

    size_t pos = fileName.find_last_of("\\/");
    if (pos != std::string::npos)
    {
      path = fileName.substr(0, pos + 1);
    }
    return path;
  }
};


template <typename Space>
struct Settings : public BasicSettings
{
  SPACE_TYPEDEFS

  MeshSettings            <Space>   mesh;
  ScheduleSettings        <Space>   schedule;
  TaskSettings            <Space>   task;
  SolverSettings          <Space>   solver;
  DetectorsSettings       <Space>   detectors;
  MeshBuilderSettings     <Space>   meshBuilder;
  ResultCombinerSettings  <Space>   resultCombiner;
  std::vector< SnapshotSettings<Space> > snapshots;

  std::string settingsFileName;
  void Parse(const char* fileName) override;

private:
  void ParseSettingsFile(const char* fileName);
};

std::string BasicSettings::ParseSettingsFileName(const char* fileName)
{
  TiXmlDocument taskFile;
  TiXmlElement* settingsElement = GetSettingXmlElement(taskFile, fileName);

  std::string settingsFileName;
  ParseString(settingsElement, "fileName", &settingsFileName);
  if (settingsFileName.empty())
  {
    std::cerr << "You should choose settings file to compute\n";
    throw;
  }
  return settingsFileName;
}

void BasicSettings::Parse(const char* fileName)
{
  std::string settingsFileName = ParseSettingsFileName(fileName);
  BasicSettings::ParseSettingsFile(settingsFileName.c_str());
}

void BasicSettings::ParseSettingsFile(const char* fileName)
{
  TiXmlDocument settingsFile;
  TiXmlElement* settingsElement = GetSettingXmlElement(settingsFile, fileName);

  if(settingsElement->QueryIntAttribute("dimsCount", &configDimsCount) != TIXML_SUCCESS)
  {
    configDimsCount = 2;
  }
}

template<typename Space>
void Settings<Space>::Parse(const char* fileName)
{
  this->settingsFileName = ParseSettingsFileName(fileName);
  Settings<Space>::ParseSettingsFile(settingsFileName.c_str());
}

template<typename Space>
void Settings<Space>::ParseSettingsFile(const char* fileName)
{
  TiXmlDocument settingsFile;
  TiXmlElement* settingsElement = GetSettingXmlElement(settingsFile, fileName);

  if(settingsElement->QueryIntAttribute("dimsCount", &configDimsCount) != TIXML_SUCCESS)
  {
    configDimsCount = 2;
  }

  TiXmlElement* meshInfoElement = settingsElement->FirstChildElement("Mesh");
  if (meshInfoElement)
  {
    mesh.Parse(meshInfoElement);
  } else
  {
    std::cout << "There is no Mesh section";
  }

  TiXmlElement* snapshotElement;
  for (IndexType snapshotIndex = 0; ; ++snapshotIndex)
  {
    if (snapshotIndex == 0)
    {
      snapshotElement = settingsElement->FirstChildElement("Snapshot");
    } else
    {
      snapshotElement = snapshotElement->NextSiblingElement("Snapshot");
    }

    SnapshotSettings<Space> snapshot;
    if (!snapshotElement)
    {
      snapshot.data      .used = false;
      snapshot.mesh      .used = false;
      snapshot.contacts  .used = false;
      if (snapshotIndex == 0)
      {
        std::cout << "There is no Snapshot section";
      }
      break;
    }else
    {
      snapshot.Parse(snapshotElement);
    }
    snapshots.push_back(snapshot);
  }

  TiXmlElement* scheduleElement = settingsElement->FirstChildElement("Schedule");
  if (scheduleElement)
  {
    schedule.Parse(scheduleElement);
  } else
  {
    std::cout << "There is no Schedule element";
  }
    
  TiXmlElement* taskElement = settingsElement->FirstChildElement("Task");
  if (taskElement)
  {
    task.Parse(taskElement);
  } else
  {
    std::cout << "There is no Task element";
  }

  TiXmlElement* solverElement = settingsElement->FirstChildElement("Solver");
  if (solverElement)
  {
    solver.Parse(solverElement);
  } else
  {
    std::cout << "There is no Solver section";
  }

  TiXmlElement* meshBuilderElement = settingsElement->FirstChildElement("MeshBuilder");
  if (meshBuilderElement)
  {
    meshBuilder.Parse(meshBuilderElement);
  } else
  {
    std::cout << "There is no MeshBuilder section";
  }

  TiXmlElement* resultCombinerElement = settingsElement->FirstChildElement("ResultCombiner");
  if (resultCombinerElement)
  {
    resultCombiner.Parse(resultCombinerElement);
  } else
  {
    std::cout << "There is no ResultCombiner section";
  }


  detectors.Parse(settingsElement);
}

#pragma once

#include "ParserUtil.h"
#include "../../../Maths/Spaces.h"

template <typename Space>
struct ResultCombinerSettings
{
  SPACE_TYPEDEFS

  ResultCombinerSettings():
    snapshotsCount(1)
  {}

  void Parse(TiXmlElement *resultCombinerElement)
  {
    snapshotsCount = 0;
    TiXmlElement* snapshotsElement = resultCombinerElement->FirstChildElement("Snapshots");
    if(snapshotsElement)
    {
      ParseUnsigned(snapshotsElement, "count", &snapshotsCount);
    }

    TiXmlElement* seismogramsElement = resultCombinerElement->FirstChildElement("Seismograms");
    if(seismogramsElement)
    {
      ParseString(seismogramsElement, "detectorsFileName", &detectorsFileName);

      ParseString(seismogramsElement, "detectorsLocationsFile", &detectorsLocationsFile);

      ParseString(seismogramsElement, "velocityCsvName", &velocityCsvName);
      ParseString(seismogramsElement, "velocityRefCsvData", &velocityRefCsvName);
      ParseString(seismogramsElement, "pressureCsvName", &pressureCsvName);

      ParseString(seismogramsElement, "velocityCoordSegyName", &velocityCoordSegyName);
      ParseString(seismogramsElement, "velocityDiffCoordSegyName", &velocityDiffCoordSegyName);
    }
  }

  int snapshotsCount;

  std::string detectorsFileName;

  std::string detectorsLocationsFile;

  std::string velocityCsvName;
  std::string velocityRefCsvName;
  std::string pressureCsvName;

  std::string velocityCoordSegyName;
  std::string velocityDiffCoordSegyName;
};
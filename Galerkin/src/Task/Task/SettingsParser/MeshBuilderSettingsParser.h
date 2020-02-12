#pragma once

#include "ParserUtil.h"
#include "../../../Maths/Spaces.h"

template <typename Space>
struct MeshBuilderSettings
{
  SPACE_TYPEDEFS

  MeshBuilderSettings(): 
    partitionAlgorithmName("METIS_PartMeshDual"),  // METIS_PartMeshNodal
    salomeFileName("none"),
    meshFileName("")
  {}

  void Parse(TiXmlElement *meshBuilderElement)
  {
    TiXmlElement* partitionInfoElement = meshBuilderElement->FirstChildElement("Partition");
    if(partitionInfoElement)
    {
      ParseString(partitionInfoElement, "algorithm", &partitionAlgorithmName);
    }
  
    TiXmlElement* salomeFileElement = meshBuilderElement->FirstChildElement("SalomeFile");

    assert(salomeFileElement);
    ParseString(salomeFileElement, "fileName", &salomeFileName);

    TiXmlElement* outFileElement = meshBuilderElement->FirstChildElement("OutFile");
    if(outFileElement)
    {
      ParseString(outFileElement, "fileName", &meshFileName);
    }
  }

  std::string partitionAlgorithmName;
  std::string salomeFileName;
  std::string meshFileName;
};
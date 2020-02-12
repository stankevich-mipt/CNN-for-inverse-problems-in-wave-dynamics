#pragma once

#include "../../Utils/Base64.h"
#include "../../3rdparty/tinyxml/tinyxml.h"
#include "../../3rdparty/tinyxml/tinystr.h"
#include <zlib.h>
#include <zconf.h>

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

using std::string;
using std::fstream;

template <typename ElasticSpace>
class VtkReader 
{
public:
  typedef typename ElasticSpace::SpaceType Space;
  SPACE_TYPEDEFS

  typedef typename ElasticSpace::Elastic   Elastic;
  
  VtkReader()
  {
    ::build_decoding_table();
  }

  ~VtkReader()
  {
    ::base64_cleanup();
  }

  // return false if error occured or format was wrong
  bool Read(const string& fileName,
            std::vector<Elastic>* const data)
  {
    TiXmlDocument file;
    bool fileLoadOkay = file.LoadFile(fileName);

    if (!fileLoadOkay)
    {
      std::cerr << "Loading " << fileName << " fails with following error: " << 
        std::string(file.ErrorDesc()) << 
        " in row " << file.ErrorRow() << std::endl;
      return false;
    }

    TiXmlElement* vtkFileElement = file.FirstChildElement("VTKFile");
    if (!vtkFileElement)
    {
      std::cerr << "There is no VTKFile element" << std::endl;
      return false;
    }

    TiXmlElement* imageDataElement = vtkFileElement->FirstChildElement("ImageData");
    if (!imageDataElement)
    {
      std::cerr << "There is no ImageData element" << std::endl;
      return false;
    }

    TiXmlElement* pieceDataElement = imageDataElement->FirstChildElement("Piece");
    if (!pieceDataElement)
    {
      std::cerr << "There is no Piece element" << std::endl;
      return false;
    }

    IndexType snapshotSize;
    if (ParseExtent(pieceDataElement, "Extent", &snapshotSize) != TIXML_SUCCESS)
    {
      std::cerr << "Can`t read Extent" << std::endl;
      return false;
    }

    TiXmlElement* pointDataElement = pieceDataElement->FirstChildElement("PointData");
    if (!pointDataElement)
    {
      std::cerr << "There is no PointData element" << std::endl;
      return false;
    }

    data->resize(snapshotSize);
    const IndexType tensionElementsCount = Space::Dimension * (Space::Dimension + 1) / 2;
    if (rawData.size() < snapshotSize * tensionElementsCount) rawData.resize(snapshotSize * tensionElementsCount);

    TiXmlElement* dataArrayElement = nullptr;
    for (IndexType dataArrayIndex = 0; dataArrayIndex < 2; ++dataArrayIndex)
    {
      if (dataArrayIndex == 0)
      {
        dataArrayElement = pointDataElement->FirstChildElement("DataArray");
      } else
      {
        dataArrayElement = dataArrayElement->NextSiblingElement("DataArray");
      }

      if (!dataArrayElement)
      {
        if (dataArrayIndex == 0)
        {
          std::cerr << "There is no DataArray element" << std::endl;
          return false;
        } else
        return true;
      }

      std::string dataName;
      if (ParseString(dataArrayElement, "Name", &dataName) != TIXML_SUCCESS)
      {
        std::cerr << "Can`t read Name" << std::endl;
        return false;
      }

      if (!ParseDataArray(dataArrayElement->GetText(), rawData))
      {
        std::cerr << "Can`t parse data array" << std::endl;
        return false;
      }

      if (dataName == "Velocity")
      {
        for (IndexType dataIndex = 0; dataIndex < snapshotSize; ++dataIndex)
        {
          SetVelocity(Overload<Space>(), dataIndex, data);
        }
      } else 
      if (dataName == "Stress")
      {
        for (IndexType dataIndex = 0; dataIndex < snapshotSize; ++dataIndex)
        {
          SetTension(Overload<Space>(), dataIndex, data);
        }
      }
    }
    return true;
  }

private:
  void SetVelocity(Overload<Space2>, IndexType dataIndex, std::vector<Elastic>* const data)
  {
    (*data)[dataIndex].SetVelocity(Vector(rawData[3 * dataIndex + 0], rawData[3 * dataIndex + 1]));
    assert(rawData[3 * dataIndex + 2] == Scalar(0.0));
  }

  void SetVelocity(Overload<Space3>, IndexType dataIndex, std::vector<Elastic>* const data)
  {
    (*data)[dataIndex].SetVelocity(Vector(rawData[3 * dataIndex + 0], rawData[3 * dataIndex + 1], rawData[3 * dataIndex + 2]));
  }

  void SetTension(Overload<Space2>, IndexType dataIndex, std::vector<Elastic>* const data)
  {
    (*data)[dataIndex].SetTension(
      rawData[3 * dataIndex + 0], 
      rawData[3 * dataIndex + 1], 
      rawData[3 * dataIndex + 2]
    );
  }

  void SetTension(Overload<Space3>, IndexType dataIndex, std::vector<Elastic>* const data)
  {
    (*data)[dataIndex].SetTension(
      rawData[6 * dataIndex + 0], 
      rawData[6 * dataIndex + 1], 
      rawData[6 * dataIndex + 2],
      rawData[6 * dataIndex + 3],
      rawData[6 * dataIndex + 4],
      rawData[6 * dataIndex + 5]
    );
  }

  int ParseExtent(TiXmlElement* element, const std::string& name, IndexType* snapshotSize)
  {
    std::string texValue;
    element->QueryStringAttribute(name.c_str(), &texValue);

    *snapshotSize = 1;
    std::stringstream stream(texValue);

    for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
    {
      int imin, imax;
      stream >> imin >> imax;
      (*snapshotSize) *= imax - imin + 1;
    }

    return stream.bad() ? TIXML_WRONG_TYPE : TIXML_SUCCESS;
  }

  int ParseString(TiXmlElement *element, const std::string name, std::string *value)
  {
    return element->QueryStringAttribute(name.c_str(), value);
  }

  bool ParseDataArray(const char* input, std::vector<Scalar>& data)
  {
    int prefixSize = 4 * ((4 * sizeof(uint32_t) + 2) / 3);

    size_t outputSize;
    unsigned char* output = ::base64_decode(input, prefixSize, &outputSize);
    if (!output)
    {
      std::cerr << "Can`t decode data array prefix" << std::endl;
      return false;
    }

    size_t sourceLen = ((uint32_t*)output)[2];
    int compressedLen = ((uint32_t*)output)[3];
    free(output);
    output = nullptr;

    int dataSize = 4 * ((compressedLen + 2) / 3);
    output = ::base64_decode(input + prefixSize, dataSize, &outputSize);
    if (!output)
    {
      std::cerr << "Can`t decode data" << std::endl;
      return false;
    }

    if (data.size() * sizeof(Scalar) < sourceLen) data.resize(sourceLen / sizeof(Scalar) + 1);
    int result = ::uncompress((Bytef*)(data.data()), (uLongf*)&sourceLen, output, compressedLen);
    
    free(output);

    return result == Z_OK;
  }

  std::vector<Scalar> rawData;
};

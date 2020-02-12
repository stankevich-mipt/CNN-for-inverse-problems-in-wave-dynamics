#pragma once

#include "../../GeomMesh/MeshIO/Local/MeshIO.h"
#include "../../../IO/Vtk/MeshVtkWriter.h"

template <typename Space>
class BasicMeshBuilder
{
public:
  SPACE_TYPEDEFS
  virtual ~BasicMeshBuilder() = default;

  virtual void BuildMesh(MeshIO<Space>* const mesh, 
    GeomMesh<Space>* geomMesh, 
    MeshIO<Space>* const nonDuplicatedVerticesMesh = nullptr) = 0;

  virtual void BuildMediumParams(MeshIO<Space>* const domain, std::vector<char>* mediumParams) {}

  void SaveGeometry(MeshIO<Space>* const mesh, const std::string& fileName)
  {
    MeshVtkWriter<Space> writer;
    writer.Write(fileName, mesh->vertices, mesh->indices);
  }

protected:
  template <typename DataType>
  void WriteData(std::vector<char>* mediumParams, DataType data)
  {
    char* bytes = static_cast<char*>(&data);
    for (IndexType byteIndex = 0; byteIndex < sizeof(DataType); ++byteIndex)
    {
      mediumParams->push_back(bytes[byteIndex]);
    }
  }
};

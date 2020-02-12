#pragma once

#include "UnstructuredMeshBuilder.h"
#include <vector>
#include <fstream>

template <typename Space>
class SalomeMeshBuilder: public UnstructuredMeshBuilder<Space>
{
public:
  SPACE_TYPEDEFS

  SalomeMeshBuilder(const std::string& meshFileName, bool duplicateNodes = false, Scalar splitNodeBiasRatio = 0.0): 
    UnstructuredMeshBuilder<Space>(duplicateNodes, splitNodeBiasRatio),
    meshFileName(meshFileName)
  {
  }

  virtual void LoadGeom(MeshIO<Space>* const mesh) override
  {
    mesh->Load(std::string(AddExtensionToFileName(meshFileName, ".mesh")).c_str(), IO::Ascii);
  }

  virtual void BuildMediumParams(MeshIO<Space>* const mesh, std::vector<char>* mediumParams) override
  {
    std::fstream paramsFile(std::string(AddExtensionToFileName(meshFileName, ".params")).c_str(), std::fstream::in | std::fstream::binary);
    
    mediumParams->resize(mesh->GetCellsCount());
    if (paramsFile.fail())
    {
      std::cerr << "Can`t read params file" << std::endl;
      std::fill(mediumParams->begin(), mediumParams->end(), 0);
      return;
    }
    paramsFile.read(mediumParams->data(), mediumParams->size());
  }

  void BuildDetectors(MeshIO<Space>* const mesh) override
  {
  }

  virtual void BuildAdditionalContacts(MeshIO<Space>* const mesh) override
  {
  }

private:
  std::string meshFileName;
};

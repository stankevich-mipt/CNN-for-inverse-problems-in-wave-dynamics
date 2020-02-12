#pragma once

#include "BasicVtkWriter.h"
#include "../../Maths/Spaces.h"
#include "../../Task/GeomMesh/GeomMesh/GeomMesh.h"
#include "../../Task/ElasticSystem/ElasticSystem.h"

#include <set>

template<typename Space, typename FunctionSpace>
class ContactVtkWriter: public BasicVtkWriter<Space>
{
public:
  SPACE_TYPEDEFS

  enum FaceType {Contact, Boundary};

  typedef typename GeomMesh<Space>::Node Node;

  using BasicVtkWriter<Space>::TypeName;

  struct OutputData
  {
    std::vector<Node> nodes;

    struct Face
    {
      Face(): faceType(0)
      {}
      IndexType incidentNodes[Space::NodesPerFace];
      IndexType faceType;
    };

    std::vector<Face> faces;
  };

  void Write(const std::string& fileName,
    ElasticVolumeMesh<Space, FunctionSpace>* mesh,
    const std::set<IndexType>& faceTypesToDraw,
    FaceType faceType) const;

protected:
  void SaveToFile(const std::string& fileName, const OutputData& outputData) const;

private:
  OutputData ConstructOutputData(ElasticVolumeMesh<Space, FunctionSpace>* mesh,
    const std::set<IndexType>& faceTypesToDraw,
    FaceType faceType) const;
};

#include "ContactVtkWriter.inl"

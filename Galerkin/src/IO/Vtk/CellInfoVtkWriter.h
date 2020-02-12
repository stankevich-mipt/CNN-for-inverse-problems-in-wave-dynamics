#include "BasicVtkWriter.h"
#include "../../Maths/Spaces.h"
#include "../../Task/GeomMesh/GeomMesh/GeomMesh.h"
#include "../../Task/ElasticSystem/ElasticSystem.h"

template<typename Space, typename FunctionSpace>
class CellInfoVtkWriter : public BasicVtkWriter<Space>
{
public:
  SPACE_TYPEDEFS

  typedef typename GeomMesh<Space>::Node              Node;
  typedef typename GeomMesh<Space>::Cell              Cell;

  using BasicVtkWriter<Space>::TypeName;

  struct OutputData
  {
    std::vector<Node> nodes;
    std::vector<Cell> cells;

    struct CellData
    {
      CellData(IndexType isCellBroken, Scalar plasticDeforms = 0) :
        isCellBroken(isCellBroken), plasticDeforms(plasticDeforms)
      {}

      CellData()
      {}

      IndexType isCellBroken;
      Scalar plasticDeforms;
      Scalar density;
    };
    std::vector<CellData> cellData;
  };

  void Write(const std::string& fileName,
    ElasticVolumeMesh<Space, FunctionSpace>* mesh) const;

protected:
  void SaveToFile(const std::string& fileName, const OutputData& outputData) const;

private:
  OutputData ConstructOutputData(ElasticVolumeMesh<Space, FunctionSpace>* mesh) const;
};

#include "CellInfoVtkWriter.inl"

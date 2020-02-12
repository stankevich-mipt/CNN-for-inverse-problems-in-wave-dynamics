#pragma once

#include "../../Maths/Spaces.h"
#include "BasicVtkWriter.h"
#include "../../Task/GeomMesh/GeomMesh/GeomMesh.h"

template < typename Space, typename CellInfo = typename AdditionalCellInfo<Space>:: template AuxInfo<char> >
class MeshVtkWriter: public BasicVtkWriter<Space>
{
  public:
    SPACE_TYPEDEFS

    typedef typename GeomMesh<Space>::Node NodeType;
    typedef typename GeomMesh<Space>::Cell CellType;

    using BasicVtkWriter<Space>::TypeName;

    void Write(const std::string& fileName,
      const std::vector<NodeType>& nodes, const std::vector<CellType>& cells,
      std::vector<typename CellInfo::DataType>& cellsData,
      std::vector<bool>* needToDrawCell = nullptr);

    void Write(const std::string& fileName,
      const std::vector<NodeType>& nodes, const std::vector<CellType>& cells,
      const std::vector<CellInfo>& cellInfos  = std::vector<CellInfo>(),
      std::vector<bool>* needToDrawCell = nullptr);

    void Write(const std::string& fileName,
      const std::vector<Vector>& nodes, 
      const std::vector<IndexType>& indices,
      std::vector<typename CellInfo::DataType>& cellsData,
      std::vector<bool>* needToDrawCell = nullptr);

    void Write(const std::string& fileName,
      const std::vector<Vector>& nodes, 
      const std::vector<IndexType>& indices,
      const std::vector<CellInfo>& cellInfos  = std::vector<CellInfo>(),
      std::vector<bool>* needToDrawCell = nullptr);

    void WriteFaces(const std::string& fileName,
      const std::vector<Vector>& nodes,
      const std::vector<IndexType>& indices,
      const std::vector<CellInfo>& cellInfos = std::vector<CellInfo>());

  private:
    void WriteCellData(std::fstream& file,
                       const std::vector<CellInfo>& cellInfos  = std::vector<CellInfo>(),
      std::vector<bool>* needToDrawCell = nullptr);

    IndexType GetCellsToDrawCount(std::vector<bool>* needToDrawCell, IndexType totalCellsCount)
    {
      IndexType cellsToDrawCount = 0;
      for (IndexType cellIndex = 0; cellIndex < totalCellsCount; ++cellIndex)
      {
        if (needToDrawCell == nullptr || (*needToDrawCell)[cellIndex]) ++cellsToDrawCount;
      }
      return cellsToDrawCount;
    }
};

template <typename Space, typename CellInfo>
void MeshVtkWriter<Space, CellInfo>::Write(const std::string& fileName,
  const std::vector<NodeType>& nodes, const std::vector<CellType>& cells,
  std::vector<typename CellInfo::DataType>& cellsData,
  std::vector<bool>* needToDrawCell)
{
  std::vector<CellInfo> cellInfos(cells.size());
  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    cellInfos[cellIndex].count = cellsData.size() > 0 ? 1 : 0;
    cellInfos[cellIndex].data = cellsData.size() > 0 ? &cellsData[cellIndex] : 0;
  }
  Write(fileName, nodes, cells, cellInfos, needToDrawCell);
}

template <typename Space, typename CellInfo>
void MeshVtkWriter<Space, CellInfo>::Write(const std::string& fileName,
  const std::vector<Vector>& nodes, 
  const std::vector<IndexType>& indices,
  std::vector<typename CellInfo::DataType>& cellsData,
  std::vector<bool>* needToDrawCell)
{
  std::vector<CellInfo> cellInfos(indices.size() / Space::NodesPerCell);
  for (IndexType cellIndex = 0; cellIndex < indices.size() / Space::NodesPerCell; ++cellIndex)
  {
    cellInfos[cellIndex].count = cellsData.size() > 0 ? 1 : 0;
    cellInfos[cellIndex].data = cellsData.size() > 0 ? &cellsData[cellIndex] : 0;
  }
  Write(fileName, nodes, indices, cellInfos, needToDrawCell);
}


template <typename Space, typename CellInfo>
void MeshVtkWriter<Space, CellInfo>::Write(const std::string& fileName,
  const std::vector<NodeType>& nodes, 
  const std::vector<CellType>& cells,
  const std::vector<CellInfo>& cellInfos,
  std::vector<bool>* needToDrawCell)
{  
  std::fstream file(fileName.c_str(), std::fstream::out | std::fstream::binary);
  
  if (file.fail()) 
  {
    return;
  }

  // header lines
  file << "# vtk DataFile Version 4.0\n";
  file << "vtk output\n";
  file << "ASCII\n";

  file << "DATASET UNSTRUCTURED_GRID\n";  

  // points
  file << "POINTS " << nodes.size() << " " << TypeName() << "\n"; 
  for (IndexType index = 0; index < nodes.size(); ++index) 
  {
    for (IndexType coordIndex = 0; coordIndex < 3; ++coordIndex)
    {
      file << nodes[index].pos.Get(coordIndex) << " ";
    }
  }
  file << "\n";

  IndexType cellsToDrawCount = GetCellsToDrawCount(needToDrawCell, cells.size());

  // cells
  file << "CELLS " << cellsToDrawCount << " "
                      << cellsToDrawCount * (Space::NodesPerCell /*vertex count in a cell*/
                                         + 1 /*count of connectivity indices*/) << "\n";
    
  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex) 
  {
    if (needToDrawCell == nullptr || (*needToDrawCell)[cellIndex])
    {
      file << Space::NodesPerCell;
      for (IndexType nodeIndex = 0; nodeIndex < Space::NodesPerCell; ++nodeIndex)
      {
        file << " " << cells[cellIndex].incidentNodes[nodeIndex];
      }
      file << "\n";
    }
  }

  file << "CELL_TYPES " << cellsToDrawCount << "\n";
  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex) 
  {
    if (needToDrawCell == nullptr || (*needToDrawCell)[cellIndex])
    {
      if (Space::NodesPerCell == 4)
      {
        file << 10; // tetrahedron
      }
      else
      if (Space::NodesPerCell == 3)
      {
        file << 5; // triangle
      }
      file << std::endl;
    }
  }

  WriteCellData(file, cellInfos, needToDrawCell);
  file.close();
} 

template <typename Space, typename CellInfo>
void MeshVtkWriter<Space, CellInfo>::Write(const std::string& fileName,
  const std::vector<Vector>& nodes, 
  const std::vector<IndexType>& indices,
  const std::vector<CellInfo>& cellInfos,
  std::vector<bool>* needToDrawCell)
{
  std::fstream file(fileName.c_str(), std::fstream::out | std::fstream::binary);
  
  if (file.fail()) 
  {
    return;
  }

  // header lines
  file << "# vtk DataFile Version 4.0\n";
  file << "vtk output\n";
  file << "ASCII\n";

  file << "DATASET UNSTRUCTURED_GRID\n";  

  // points
  file << "POINTS " << nodes.size() << " " << TypeName() << "\n"; 
  for (IndexType index = 0; index < nodes.size(); ++index) 
  {
    for (IndexType coordIndex = 0; coordIndex < 3; ++coordIndex)
    {
      file << nodes[index].Get(coordIndex) << " ";
    }
  }
  file << "\n";

  IndexType cellsToDrawCount = GetCellsToDrawCount(needToDrawCell, indices.size() / Space::NodesPerCell);

  // cells
  file << "CELLS " << cellsToDrawCount << " "
                      << cellsToDrawCount * (Space::NodesPerCell /*vertex count in a cell*/
                                                                 + 1 /*count of connectivity indices*/) << "\n";
    
  for (IndexType cellIndex = 0; Space::NodesPerCell * cellIndex < indices.size(); ++cellIndex) 
  {
    if (needToDrawCell == nullptr || (*needToDrawCell)[cellIndex])
    {
      file << Space::NodesPerCell;
      for (IndexType nodeIndex = 0; nodeIndex < Space::NodesPerCell; ++nodeIndex)
      {
        file << " " << indices[Space::NodesPerCell * cellIndex + nodeIndex];
      }
      file << "\n";
    }
  }

  file << "CELL_TYPES " << cellsToDrawCount << "\n";
  for (IndexType cellIndex = 0; Space::NodesPerCell * cellIndex < indices.size(); ++cellIndex)
  {
    if (needToDrawCell == nullptr || (*needToDrawCell)[cellIndex])
    {
      if (Space::NodesPerCell == 4)
      {
        file << 10; // tetrahedron
      }
      else
      if (Space::NodesPerCell == 3)
      {
        file << 5; // triangle
      }
      file << std::endl;
    }
  }

  WriteCellData(file, cellInfos, needToDrawCell);
  file.close();
}

template <typename Space, typename CellInfo>
void MeshVtkWriter<Space, CellInfo>::WriteFaces(const std::string& fileName,
  const std::vector<Vector>& nodes, const std::vector<IndexType>& indices, const std::vector<CellInfo>& cellInfos)
{
  std::fstream file(fileName.c_str(), std::fstream::out | std::fstream::binary);

  if (file.fail())
  {
    return;
  }

  // header lines
  file << "# vtk DataFile Version 4.0\n";
  file << "vtk output\n";
  file << "ASCII\n";

  file << "DATASET UNSTRUCTURED_GRID\n";

  // points
  file << "POINTS " << nodes.size() << " " << TypeName() << "\n";
  for (IndexType index = 0; index < nodes.size(); ++index)
  {
    for (IndexType coordIndex = 0; coordIndex < 3; ++coordIndex)
    {
      file << nodes[index].Get(coordIndex) << " ";
    }
  }
  file << "\n";

  // cells
  file << "CELLS " << indices.size() / Space::NodesPerFace << " "
    << indices.size() / Space::NodesPerFace * (Space::NodesPerFace /*vertex count in a cell*/
    + 1 /*count of connectivity indices*/) << "\n";

  for (IndexType faceIndex = 0; Space::NodesPerFace * faceIndex < indices.size(); ++faceIndex) {
    file << Space::NodesPerFace;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber) {
      file << " " << indices[Space::NodesPerFace * faceIndex + nodeNumber];
    }
    file << "\n";
  }

  file << "CELL_TYPES " << indices.size() / Space::NodesPerFace << "\n";
  for (IndexType cellIndex = 0; cellIndex < indices.size() / Space::NodesPerFace; ++cellIndex)
  {
    if (Space::NodesPerFace == 3)
    {
      file << 5; // triangle
    } else
    if (Space::NodesPerFace == 2)
    {
      file << 3; // line
    }
    file << std::endl;
  }

  WriteCellData(file, cellInfos);
  file.close();
}

template <typename Space, typename CellInfo>
void MeshVtkWriter<Space, CellInfo>::WriteCellData(std::fstream& file,
  const std::vector<CellInfo>& cellInfos,
  std::vector<bool>* needToDrawCell)
{
  // cells data
  IndexType dataSize = 0;
  if (!cellInfos.empty())
  {
    dataSize = cellInfos[0].count;
    for (IndexType index = 1; index < cellInfos.size(); ++index)
    {
      dataSize = std::max(dataSize, cellInfos[index].count);
    }
  }

  if (dataSize && dataSize != IndexType(-1))
  {
    file << std::endl;

    IndexType cellsToDrawCount = GetCellsToDrawCount(needToDrawCell, cellInfos.size());

    file << "CELL_DATA " << cellsToDrawCount << std::endl;
    file << "SCALARS cell_scalars " << TypeName() << " " << dataSize << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;

    for (IndexType cellIndex = 0; cellIndex < cellInfos.size(); ++cellIndex)
    {
      if (needToDrawCell == nullptr || (*needToDrawCell)[cellIndex])
      {
        for (IndexType dataIndex = 0; dataIndex < dataSize; ++dataIndex)
        {
          file << (cellInfos[cellIndex].count > dataIndex ? cellInfos[cellIndex].data[dataIndex] : 0) << " ";
        }
        file << std::endl;
      }
    }
  }
}

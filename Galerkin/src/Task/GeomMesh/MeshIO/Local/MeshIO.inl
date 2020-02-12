#include <algorithm>

template <typename Space>
MeshIO<Space>::MeshIO(): 
  contactTypesCount(0),
  boundaryTypesCount(0)
{
}

template <typename Space>
void MeshIO<Space>::Save(const std::string& fileName, IO::FileType fileType)
{
  std::fstream file;

  switch (fileType)
  {
    case IO::Ascii: file.open(fileName.c_str(), std::fstream::out); break;
    case IO::Binary: file.open(fileName.c_str(), std::fstream::out | std::fstream::binary); break;
  }
  
  if (file.fail())
  {
    std::cerr << "Can`t open file " << fileName << " for saving" << std::endl;
    return;
  }

  Save(file, fileType);
  file.close();
}

template <typename Space>
void MeshIO<Space>::Load(const std::string& fileName, IO::FileType fileType)
{
  std::fstream file;

  switch (fileType)
  {
    case IO::Ascii: file.open(fileName.c_str(), std::fstream::in); break;
    case IO::Binary: file.open(fileName.c_str(), std::fstream::in | std::fstream::binary); break;
  }
  
  if (file.fail())
  {
    std::cerr << "Can`t open file " << fileName << " for loading" << std::endl;
    throw;
  }

  Load(file, fileType);
  file.close();
}

template <typename Space>
typename Space::IndexType MeshIO<Space>::GetCellsCount() const
{
  return indices.size() / Space::NodesPerCell;
}

template <typename Space>
typename Space::Vector MeshIO<Space>::GetCellCenter(IndexType cellIndex) const
{
  Vector center = Vector::zero();
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    center += vertices[indices[cellIndex * Space::NodesPerCell + nodeNumber]];
  }
  center /= Scalar(Space::NodesPerCell);
  return center;
}

template <typename Space>
bool MeshIO<Space>::operator==(const MeshIO<Space>& other) const
{
  if (vertices.size() != other.vertices.size()) return false;

  const Scalar eps = 1e-3;
  for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
  {
    Scalar meanLen = (vertices[vertexIndex] + other.vertices[vertexIndex]).Len() * Scalar(0.5);
    if ((vertices[vertexIndex] - other.vertices[vertexIndex]).Len() > meanLen * eps) return false;
  }

  if (detectorsPositions.size() != other.detectorsPositions.size()) return false;
  for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
  {
    Scalar meanLen = (detectorsPositions[detectorIndex] + other.detectorsPositions[detectorIndex]).Len() * Scalar(0.5);
    if ((detectorsPositions[detectorIndex] - other.detectorsPositions[detectorIndex]).Len() > meanLen * eps) return false;
  }

  return MeshIOBase<Space>::operator==(other) && 
    indices == other.indices &&
    contactTypesCount == other.contactTypesCount &&
    boundaryTypesCount == other.boundaryTypesCount;
}

// file should be opened
template <typename Space>
void MeshIO<Space>::Save(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      // vertices
      file << vertices.size() << std::endl;
      for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file << vertices[vertexIndex][componentNumber] << " ";
        }
        file << std::endl;
      }
      file << std::endl;

      // indices
      file << indices.size() << std::endl;
      for (IndexType index = 0; index < indices.size(); ++index)
      {
        file << indices[index] << " ";
      }
      file << std::endl;

      // contacts
      SaveContacts(file, IO::Ascii);
      std::cout << std::endl;

      // boundaries
      SaveBoundaries(file, IO::Ascii);
      std::cout << std::endl;

      // detectors
      file << detectorsPositions.size() << std::endl;
      for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file << detectorsPositions[detectorIndex][componentNumber] << " ";
        }
        file << std::endl;
      }
      file << std::endl;
    break;
    case IO::Binary:
      // vertices
      IO::Write(file, vertices.size());
      IO::WriteVector(file, vertices);

      IO::Write(file, indices.size());
      IO::WriteVector(file, indices);

      // contacts
      SaveContacts(file, IO::Binary);

      // boundaries
      SaveBoundaries(file, IO::Binary);

      // detectors
      IO::Write(file, detectorsPositions.size());
      IO::WriteVector(file, detectorsPositions);
    break;
  }
}

template <typename Space>
void MeshIO<Space>::Load(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
    {
      // vertices
      IndexType verticesCount;
      file >> verticesCount;
      vertices.resize(verticesCount);

      for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file >> vertices[vertexIndex][componentNumber];
        }
      }

      // indices
      IndexType indicesCount;
      file >> indicesCount;
      indices.resize(indicesCount);
      for (IndexType index = 0; index < indices.size(); ++index)
      {
        file >> indices[index];
      }

      // contacts
      LoadContacts(file, IO::Ascii);

      // boundaries
      LoadBoundaries(file, IO::Ascii);

      // detectors
      IndexType detectorsCount;
      file >> detectorsCount;
      detectorsPositions.resize(detectorsCount);
      for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file >> detectorsPositions[detectorIndex][componentNumber];
        }
      }
    } break;

    case IO::Binary:
    {
      // vertices
      IndexType verticesCount;
      IO::Read(file, verticesCount);
      vertices.resize(verticesCount);
      IO::Read(file, vertices.data(), verticesCount);

      // indices
      IndexType indicesCount;
      IO::Read(file, indicesCount);
      indices.resize(indicesCount);
      IO::Read(file, indices.data(), indicesCount);

      // contacts
      LoadContacts(file, IO::Binary);

      // boundaries
      LoadBoundaries(file, IO::Binary);

      // detectors
      IndexType detectorsCount;
      IO::Read(file, detectorsCount);
      detectorsPositions.resize(detectorsCount);
      IO::Read(file, detectorsPositions.data(), detectorsCount);
    } break;
  }
}

template <typename Space>
void MeshIO<Space>::SaveDetectors(const std::string& vtkFileName)
{
  Scalar minY = vertices[0].y;
  Scalar maxY = vertices[0].y;
  for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
  {
    minY = std::min(minY, vertices[vertexIndex].y);
    maxY = std::max(maxY, vertices[vertexIndex].y);
  }

  std::vector<Vector> vertices;
  std::vector<IndexType> indices;

  const Scalar Height = (maxY - minY) / 50;

  for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
  {
    vertices.push_back(detectorsPositions[detectorIndex]);
    vertices.push_back(detectorsPositions[detectorIndex] + Vector::yAxis() * Height);
    indices.push_back(2 * detectorIndex + 0);
    indices.push_back(2 * detectorIndex + 1);
  }

  MeshVtkWriter<Space> writer;
  writer.WriteFaces(vtkFileName, vertices, indices);
}


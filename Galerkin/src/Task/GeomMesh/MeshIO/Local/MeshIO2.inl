template <>
void MeshIO<Space2>::SaveContacts(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      using std::endl;
      file << contactEdges.size() << endl;
      for (IndexType contactEdgesIndex = 0; contactEdgesIndex < contactEdges.size(); ++contactEdgesIndex)
      {
        for (IndexType edgePairNumber = 0; edgePairNumber < 2; ++edgePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
          {
            file << contactEdges[contactEdgesIndex].edges[edgePairNumber].nodeIndices[nodeNumber] << " ";
          }
        }
      }
      file << endl;

      assert(contactEdgesCount.size() == contactTypesCount);

      file << contactEdgesCount.size() << endl;
      for (IndexType index = 0; index < contactEdgesCount.size(); ++index)
      {
        file << contactEdgesCount[index] << " ";
      } 
      file << endl;
    break;
    case IO::Binary:
      IO::Write(file, contactEdges.size());
      IO::WriteVector(file, contactEdges);

      assert(contactEdgesCount.size() == contactTypesCount);
      IO::Write(file, contactEdgesCount.size());
      IO::WriteVector(file, contactEdgesCount);
    break;
  }
}

template <>
void MeshIO<Space2>::SaveBoundaries(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      using std::endl;
      file << boundaryEdges.size() << endl;
      for (IndexType boundaryEdgesIndex = 0; boundaryEdgesIndex < boundaryEdges.size(); ++boundaryEdgesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
        {
          file << boundaryEdges[boundaryEdgesIndex].nodeIndices[nodeNumber] << " ";
        }
      }
      file << endl;

      assert(boundaryEdgesCount.size() == boundaryTypesCount);
      file << boundaryEdgesCount.size() << endl;
      for (IndexType index = 0; index < boundaryEdgesCount.size(); ++index)
      {
        file << boundaryEdgesCount[index] << " ";
      } 
      file << endl;
    break;
    case IO::Binary:
      IO::Write(file, boundaryEdges.size());
      IO::WriteVector(file, boundaryEdges);

      assert(boundaryEdgesCount.size() == boundaryTypesCount);
      IO::Write(file, boundaryEdgesCount.size());
      IO::WriteVector(file, boundaryEdgesCount);
    break;
  }
}

template <>
void MeshIO<Space2>::LoadContacts(std::fstream& file, IO::FileType fileType)
{
  IndexType contactEdgesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> contactEdgesSize;
      contactEdges.resize(contactEdgesSize);

      for (IndexType contactEdgesIndex = 0; contactEdgesIndex < contactEdges.size(); ++contactEdgesIndex)
      {
        for (IndexType edgePairNumber = 0; edgePairNumber < 2; ++edgePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
          {
            file >> contactEdges[contactEdgesIndex].edges[edgePairNumber].nodeIndices[nodeNumber];
          }
        }
      }

      file >> contactTypesCount;
      contactEdgesCount.resize(contactTypesCount);
      for (IndexType index = 0; index < contactEdgesCount.size(); ++index)
      {
        file >> contactEdgesCount[index];
      }
    break;
    case IO::Binary:
      IO::Read(file, contactEdgesSize);
      contactEdges.resize(contactEdgesSize);
      IO::Read(file, contactEdges.data(), contactEdgesSize);

      IO::Read(file, contactTypesCount);
      contactEdgesCount.resize(contactTypesCount);
      IO::Read(file, contactEdgesCount.data(), contactEdgesCount.size());
    break;
  }
}

template <>
void MeshIO<Space2>::LoadBoundaries(std::fstream& file, IO::FileType fileType)
{
  IndexType boundaryEdgesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> boundaryEdgesSize;
      boundaryEdges.resize(boundaryEdgesSize);

      for (IndexType boundaryEdgesIndex = 0; boundaryEdgesIndex < boundaryEdges.size(); ++boundaryEdgesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
        {
          file >> boundaryEdges[boundaryEdgesIndex].nodeIndices[nodeNumber];
        }
      }

      file >> boundaryTypesCount;
      boundaryEdgesCount.resize(boundaryTypesCount);
      for (IndexType index = 0; index < boundaryEdgesCount.size(); ++index)
      {
        file >> boundaryEdgesCount[index];
      }
    break;
    case IO::Binary:
      IO::Read(file, boundaryEdgesSize);
      boundaryEdges.resize(boundaryEdgesSize);
      IO::Read(file, boundaryEdges.data(), boundaryEdgesSize);

      IO::Read(file, boundaryTypesCount);
      boundaryEdgesCount.resize(boundaryTypesCount);
      IO::Read(file, boundaryEdgesCount.data(), boundaryTypesCount);
    break;
  }
}

template <>
void MeshIO<Space2>::SaveContacts(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space2>::AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> contactTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType contactType = 0; contactType < contactTypesCount; ++contactType)
  {
    for (IndexType contactNumber = 0; contactNumber < contactEdgesCount[contactType]; ++contactNumber)
    {
      EdgePairIndices contactEdge = contactEdges[offset + contactNumber];
      for (IndexType edgeNumber = 0; edgeNumber < 2; ++edgeNumber)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
        {
          indices.push_back(contactEdge.edges[edgeNumber].nodeIndices[nodeNumber]);
        }
        contactTypes.push_back(contactType);
      }
    }
    offset += contactEdgesCount[contactType];
  }
  cellInfos.resize(contactTypes.size());
  for (IndexType contactIndex = 0; contactIndex < cellInfos.size(); ++contactIndex)
  {
    cellInfos[contactIndex].count = 1;
    cellInfos[contactIndex].data = &(contactTypes[contactIndex]);
  }
  MeshVtkWriter<Space2, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}

template <>
void MeshIO<Space2>::SaveBoundaries(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space2>::AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> boundaryTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType boundaryType = 0; boundaryType < boundaryTypesCount; ++boundaryType)
  {
    for (IndexType boundaryNumber = 0; boundaryNumber < boundaryEdgesCount[boundaryType]; ++boundaryNumber)
    {
      BoundaryEdge boundaryEdge = boundaryEdges[offset + boundaryNumber];
      for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
      {
        indices.push_back(boundaryEdge.nodeIndices[nodeNumber]);
      }
      boundaryTypes.push_back(boundaryType);
    }
    offset += boundaryEdgesCount[boundaryType];
  }
  cellInfos.resize(boundaryTypes.size());
  for (IndexType boundaryIndex = 0; boundaryIndex < cellInfos.size(); ++boundaryIndex)
  {
    cellInfos[boundaryIndex].count = 1;
    cellInfos[boundaryIndex].data = &(boundaryTypes[boundaryIndex]);
  }
  MeshVtkWriter<Space2, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}

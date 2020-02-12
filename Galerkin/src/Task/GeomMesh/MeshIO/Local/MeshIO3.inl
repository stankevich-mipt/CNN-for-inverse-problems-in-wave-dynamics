template <>
void MeshIO<Space3>::SaveContacts(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      file << contactFaces.size() << std::endl;
      for (IndexType contactFacesIndex = 0; contactFacesIndex < contactFaces.size(); ++contactFacesIndex)
      {
        for (IndexType facePairNumber = 0; facePairNumber < 2; ++facePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
          {
            file << contactFaces[contactFacesIndex].faces[facePairNumber].nodeIndices[nodeNumber] << " ";
          }
        }
      }
      file << std::endl;

      assert(contactFacesCount.size() == contactTypesCount);

      file << contactFacesCount.size() << std::endl;
      for (IndexType contactTypeIndex = 0; contactTypeIndex < contactFacesCount.size(); ++contactTypeIndex)
      {
        file << contactFacesCount[contactTypeIndex] << " ";
      } 
      file << std::endl;
    break;
    case IO::Binary:
      IO::Write(file, contactFaces.size());
      IO::WriteVector(file, contactFaces);
      
      assert(contactTypesCount == contactFacesCount.size());
      IO::Write(file, contactTypesCount);
      IO::WriteVector(file, contactFacesCount);
    break;
  }
}

template <>
void MeshIO<Space3>::SaveBoundaries(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      file << boundaryFaces.size() << std::endl;
      for (IndexType boundaryFacesIndex = 0; boundaryFacesIndex < boundaryFaces.size(); ++boundaryFacesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
        {
          file << boundaryFaces[boundaryFacesIndex].nodeIndices[nodeNumber] << " ";
        }
      }
      file << std::endl;

      assert(boundaryFacesCount.size() == boundaryTypesCount);
      file << boundaryFacesCount.size() << std::endl;
      for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryFacesCount.size(); ++boundaryTypeIndex)
      {
        file << boundaryFacesCount[boundaryTypeIndex] << " ";
      } 
      file << std::endl;
    break;
    case IO::Binary:
      IO::Write(file, boundaryFaces.size());
      IO::WriteVector(file, boundaryFaces);

      assert(boundaryFacesCount.size() == boundaryTypesCount);
      IO::Write(file, boundaryFacesCount.size());
      IO::WriteVector(file, boundaryFacesCount);
    break;
  }
}

template <>
void MeshIO<Space3>::LoadContacts(std::fstream& file, IO::FileType fileType)
{
  IndexType contactFacesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> contactFacesSize;
      contactFaces.resize(contactFacesSize);

      for (IndexType contactFacesIndex = 0; contactFacesIndex < contactFaces.size(); ++contactFacesIndex)
      {
        for (IndexType facePairNumber = 0; facePairNumber < 2; ++facePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
          {
            file >> contactFaces[contactFacesIndex].faces[facePairNumber].nodeIndices[nodeNumber];
          }
        }
      }

      file >> contactTypesCount;
      contactFacesCount.resize(contactTypesCount);
      for (IndexType contactTypeIndex = 0; contactTypeIndex < contactFacesCount.size(); ++contactTypeIndex)
      {
        file >> contactFacesCount[contactTypeIndex];
      }
    break;
    case IO::Binary:
      IO::Read(file, contactFacesSize);
      contactFaces.resize(contactFacesSize);
      IO::Read(file, contactFaces.data(), contactFacesSize);

      IO::Read(file, contactTypesCount);
      contactFacesCount.resize(contactTypesCount);
      IO::Read(file, contactFacesCount.data(), contactFacesCount.size());
    break;
  }
}

template <>
void MeshIO<Space3>::LoadBoundaries(std::fstream& file, IO::FileType fileType)
{
  IndexType boundaryFacesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> boundaryFacesSize;
      boundaryFaces.resize(boundaryFacesSize);

      for (IndexType boundaryFacesIndex = 0; boundaryFacesIndex < boundaryFaces.size(); ++boundaryFacesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
        {
          file >> boundaryFaces[boundaryFacesIndex].nodeIndices[nodeNumber];
        }
      }

      file >> boundaryTypesCount;
      boundaryFacesCount.resize(boundaryTypesCount);
      for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryFacesCount.size(); ++boundaryTypeIndex)
      {
        file >> boundaryFacesCount[boundaryTypeIndex];
      }
    break;
    case IO::Binary:
      IO::Read(file, boundaryFacesSize);
      boundaryFaces.resize(boundaryFacesSize);
      IO::Read(file, boundaryFaces.data(), boundaryFacesSize);

      IO::Read(file, boundaryTypesCount);
      boundaryFacesCount.resize(boundaryTypesCount);
      IO::Read(file, boundaryFacesCount.data(), boundaryTypesCount);
    break;
  }
}

template <>
void MeshIO<Space3>::SaveContacts(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space3>::AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> contactTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType contactType = 0; contactType < contactTypesCount; ++contactType)
  {
    for (IndexType contactNumber = 0; contactNumber < contactFacesCount[contactType]; ++contactNumber)
    {
      FacePairIndices contactFace = contactFaces[offset + contactNumber];
      for (IndexType faceNumber = 0; faceNumber < 2; ++faceNumber)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
        {
          indices.push_back(contactFace.faces[faceNumber].nodeIndices[nodeNumber]);
        }
        contactTypes.push_back(contactType);
      }
    }
    offset += contactFacesCount[contactType];
  }
  cellInfos.resize(contactTypes.size());
  for (IndexType contactIndex = 0; contactIndex < cellInfos.size(); ++contactIndex)
  {
    cellInfos[contactIndex].count = 1;
    cellInfos[contactIndex].data = &(contactTypes[contactIndex]);
  }
  MeshVtkWriter<Space3, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}

template <>
void MeshIO<Space3>::SaveBoundaries(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space3>:: AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> boundaryTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType boundaryType = 0; boundaryType < boundaryTypesCount; ++boundaryType)
  {
    for (IndexType boundaryNumber = 0; boundaryNumber < boundaryFacesCount[boundaryType]; ++boundaryNumber)
    {
      BoundaryFace boundaryFace = boundaryFaces[offset + boundaryNumber];
      for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
      {
        indices.push_back(boundaryFace.nodeIndices[nodeNumber]);
      }
      boundaryTypes.push_back(boundaryType);
    }
    offset += boundaryFacesCount[boundaryType];
  }
  cellInfos.resize(boundaryTypes.size());
  for (IndexType boundaryIndex = 0; boundaryIndex < cellInfos.size(); ++boundaryIndex)
  {
    cellInfos[boundaryIndex].count = 1;
    cellInfos[boundaryIndex].data = &(boundaryTypes[boundaryIndex]);
  }
  MeshVtkWriter<Space3, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}

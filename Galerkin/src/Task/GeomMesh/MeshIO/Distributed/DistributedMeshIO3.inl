void DistributedMeshIO<Space3>::SaveTransitionInfo(std::fstream& file, IO::FileType fileType, const std::vector< TransitionInfo >& info)
{  
  for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
  {
    DistributedMeshIOCommon<Space3>::Save<IndexType>(file, fileType, info[domainIndex].cellsIndices);
    DistributedMeshIOCommon<Space3>::Save<IndexType>(file, fileType, info[domainIndex].nodesIndices);
    this->template SaveFaces<FaceIndices>(file, fileType, info[domainIndex].facesIndices);
  }
}

void DistributedMeshIO<Space3>::LoadTransitionInfo(std::fstream& file, IO::FileType fileType, std::vector< TransitionInfo >* const info)
{
  for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
  {
    DistributedMeshIOCommon::Load<IndexType>(file, fileType, &((*info)[domainIndex].cellsIndices));
    DistributedMeshIOCommon::Load<IndexType>(file, fileType, &((*info)[domainIndex].nodesIndices));
    this->template LoadFaces<FaceIndices>(file, fileType, &((*info)[domainIndex].facesIndices));
  }
}

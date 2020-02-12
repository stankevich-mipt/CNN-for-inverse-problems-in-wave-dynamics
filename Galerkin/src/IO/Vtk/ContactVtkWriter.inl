template <typename Space, typename FunctionSpace>
void ContactVtkWriter<Space, FunctionSpace>::Write(const std::string& fileName,
  ElasticVolumeMesh<Space, FunctionSpace>* mesh,
  const std::set<IndexType>& faceTypesToDraw,
  FaceType faceType) const
{
  OutputData outputData = ConstructOutputData(mesh, faceTypesToDraw, faceType);
  SaveToFile(fileName, outputData);
}

template<typename Space, typename FunctionSpace>
typename ContactVtkWriter<Space, FunctionSpace>::OutputData
ContactVtkWriter<Space, FunctionSpace>::ConstructOutputData(
  ElasticVolumeMesh<Space, FunctionSpace>* mesh,
  const std::set<IndexType>& faceTypesToDraw,
  FaceType faceType) const
{
  OutputData outputData;

  for (IndexType nodeIndex = 0; nodeIndex < mesh->volumeMesh.nodes.size(); ++nodeIndex)
  {
    Node node(mesh->volumeMesh.nodes[nodeIndex].pos * mesh->velocityDimensionlessMult);
    outputData.nodes.push_back(node);
  }

  for (IndexType cellIndex = 0; cellIndex < mesh->volumeMesh.cells.size(); ++cellIndex)
  {
    if (mesh->volumeMesh.isCellAvailable[cellIndex])
    {
      for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
      {
        IndexType correspondingCellIndex = mesh->volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);
        IndexType correspondingFaceNumber = mesh->volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber);
        IndexType interactionType = mesh->volumeMesh.GetInteractionType(cellIndex, faceNumber);

        bool boundaryToDraw = correspondingCellIndex == IndexType(-1) && correspondingFaceNumber == IndexType(-1) && interactionType != IndexType(-1) &&
          faceTypesToDraw.find(interactionType) != faceTypesToDraw.end();

        bool contactToDraw = correspondingCellIndex != IndexType(-1) && correspondingFaceNumber != IndexType(-1) && interactionType != IndexType(-1) &&
          faceTypesToDraw.find(interactionType) != faceTypesToDraw.end();

        if ((boundaryToDraw && faceType == Boundary) || (contactToDraw && faceType == Contact))
        {
          IndexType faceNodeNumbers[Space::NodesPerFace];
          mesh->volumeMesh.GetCellFaceNodes(cellIndex, faceNumber, faceNodeNumbers);

          typename OutputData::Face face;
          for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
          {
            face.incidentNodes[nodeNumber] = faceNodeNumbers[nodeNumber];
          }
          face.faceType = interactionType;

          outputData.faces.push_back(face);
        }
      }
    }
  }
  return outputData;
}

template <typename Space, typename FunctionSpace>
void ContactVtkWriter<Space, FunctionSpace>::SaveToFile(const std::string& fileName, const OutputData& outputData) const
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
  file << "POINTS " << outputData.nodes.size() << " " << TypeName() << "\n";
  for (IndexType index = 0; index < outputData.nodes.size(); ++index)
  {
    for (IndexType coordIndex = 0; coordIndex < 3; ++coordIndex)
    {
      file << outputData.nodes[index].pos.Get(coordIndex) << " ";
    }
  }
  file << "\n";

  // cells
  file << "CELLS " << outputData.faces.size() << " "
    << outputData.faces.size() * (Space::NodesPerFace /*vertex count in a cell*/
      + 1 /*count of connectivity indices*/) << "\n";

  for (IndexType faceIndex = 0; faceIndex < outputData.faces.size(); ++faceIndex)
  {
    file << Space::NodesPerFace;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
    {
      file << " " << outputData.faces[faceIndex].incidentNodes[nodeNumber];
    }
    file << "\n";
  }

  file << "CELL_TYPES " << outputData.faces.size() << "\n";
  for (IndexType cellIndex = 0; cellIndex < outputData.faces.size(); ++cellIndex)
  {
    if (Space::NodesPerFace == 3)
    {
      file << 5; // face
    } else
    if (Space::NodesPerFace == 2)
    {
      file << 3; // edge
    }
    file << std::endl;
  }

  file << "CELL_DATA " << outputData.faces.size() << "\n";
  file << "SCALARS TYPE INT 1\n";
  file << "LOOKUP_TABLE default\n";

  for (IndexType faceIndex = 0; faceIndex < outputData.faces.size(); ++faceIndex)
  {
    file << outputData.faces[faceIndex].faceType << std::endl;
  }

  file.close();
}

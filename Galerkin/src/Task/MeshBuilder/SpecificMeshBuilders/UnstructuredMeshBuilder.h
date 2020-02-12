#pragma once

#include "../../GeomMesh/GeomMesh/GeomMesh.h"
#include "BasicMeshBuilder.h"

#include <map>
#include <vector>
#include <string.h>

template <typename Space>
class UnstructuredMeshBuilder;

template <typename Space>
class UnstructuredMeshBuilderCommon: public BasicMeshBuilder<Space>
{
public:
  SPACE_TYPEDEFS

  using BasicMeshBuilder<Space>::WriteData;

  UnstructuredMeshBuilderCommon(bool duplicateNodes, Scalar splitNodeBiasRatio): 
    duplicateNodes(duplicateNodes), splitNodeBiasRatio(splitNodeBiasRatio), minEdgeLen(std::numeric_limits<Scalar>::max())
  {
  }

  virtual ~UnstructuredMeshBuilderCommon() = default;

  bool duplicateNodes;
  Scalar splitNodeBiasRatio;
  Scalar minEdgeLen;

protected:
  GeomMesh<Space>* updatedGeomMesh;
  GeomMesh<Space>* geomMesh;
  std::map<IndexType, bool> usedCells;
  std::vector<IndexType> originalNodeIndices;

  void BuildMesh(MeshIO<Space>* const mesh, 
    GeomMesh<Space>* updatedGeomMesh, 
    MeshIO<Space>* const nonDuplicatedVerticesMesh = nullptr) override
  {
    this->updatedGeomMesh = updatedGeomMesh;
    printf("  Building triangulation...\n");
    LoadGeom(mesh);

    if (nonDuplicatedVerticesMesh)
    {
      nonDuplicatedVerticesMesh->vertices = mesh->vertices;
      nonDuplicatedVerticesMesh->indices  = mesh->indices;
    }

    printf("  Adding cavities...\n");
    AddCavities(mesh);
    printf("  Building contacts...\n");
    BuildContacts(mesh);
    printf("  Building boundaries...\n");
    BuildBoundary(mesh);

    // all unnecessary allocated memory should be free immediately
    delete geomMesh;

    std::cout << "  ";
    // this method should be call only after building contacts and boundaries
    BuildAdditionalTopology(updatedGeomMesh, mesh);

    printf("  Building detectors...\n");
    BuildDetectors(mesh);
    printf("  Building additional contacts...\n");
    BuildAdditionalContacts(mesh);
    printf("\n\n");
    updatedGeomMesh = this->updatedGeomMesh;
  }

  virtual void LoadGeom(MeshIO<Space>* const) = 0;
  virtual void BuildBoundary(MeshIO<Space>* const mesh) = 0;
  virtual void BuildContacts(MeshIO<Space>* const mesh) = 0;

  virtual void BuildDetectors(MeshIO<Space>* const)
  {
    // without detectors
  }
  virtual void BuildAdditionalContacts(MeshIO<Space>* const)
  {
    // TODO: it`s about cycled meshes
  }

  virtual void BuildCellMediumParameters(MeshIO<Space>* const mesh)
  {
  }
private:
  void AddCavities(MeshIO<Space>* const mesh)
  {
    geomMesh = new GeomMesh<Space>();
    printf("  ");
    geomMesh->LoadGeom(mesh->vertices.data(), mesh->indices.data(),
      mesh->vertices.size(), mesh->GetCellsCount());
    geomMesh->BuildTopologyInfos();
    BuildAdditionalTopology(geomMesh, mesh);

    geomMesh->BuildOriginalCellIndices(mesh->indices.data(), mesh->GetCellsCount());

    for (IndexType edgeIndex = 0; edgeIndex < geomMesh->edges.size(); ++edgeIndex)
    {
      const typename GeomMesh<Space>::Edge& edge = geomMesh->edges[edgeIndex];
      Vector edgeNodes[2];
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
      {
        edgeNodes[nodeNumber] = geomMesh->nodes[edge.incidentNodes[nodeNumber]].pos;
      }
      minEdgeLen = std::min(minEdgeLen, (edgeNodes[0] - edgeNodes[1]).Len());
    }

    for (IndexType nodeIndex = 0; nodeIndex < geomMesh->nodes.size(); ++nodeIndex)
    {
      originalNodeIndices.push_back(nodeIndex);
    }

    for (IndexType nodeIndex = 0; nodeIndex < geomMesh->nodes.size(); ++nodeIndex)
    {
      std::map<IndexType, IndexType> incidentCellsMarkers;
      IndexType cellMarkersCount;
      GetIncidentCellsMarkers(nodeIndex, incidentCellsMarkers, cellMarkersCount);
      std::vector<Vector> nodesDelta(cellMarkersCount, Vector::zero());
      // splitting node
      for (IndexType incidentCellNumber = 0; incidentCellNumber < geomMesh->nodeTopologyInfos[nodeIndex].incidentCellsCount; ++incidentCellNumber)
      {
        IndexType updatedCellIndex = geomMesh->nodeTopologyInfos[nodeIndex].incidentCells[incidentCellNumber];
        IndexType cellMarker = incidentCellsMarkers[updatedCellIndex];
        IndexType incidentCellIndex = geomMesh->originalCellIndices[updatedCellIndex];
        if (cellMarkersCount > 1)
        {
          nodesDelta[cellMarker] += (mesh->GetCellCenter(incidentCellIndex) - mesh->vertices[nodeIndex]).GetNorm();
        }
        if (cellMarker > 0)
        {
          bool nodeFound = false;
          for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
          {
            if (mesh->indices[incidentCellIndex * Space::NodesPerCell + nodeNumber] == nodeIndex)
            {
              mesh->indices[incidentCellIndex * Space::NodesPerCell + nodeNumber] = mesh->vertices.size() + cellMarker - 1;
              nodeFound = true;
              break;
            }
          }
          assert(nodeFound);
        }
      }

      mesh->vertices[nodeIndex] += nodesDelta[0].GetNorm() * minEdgeLen * splitNodeBiasRatio;
      for (IndexType duplicatedNodeNumber = 1; duplicatedNodeNumber < cellMarkersCount; ++duplicatedNodeNumber)
      {
        mesh->vertices.push_back(geomMesh->nodes[nodeIndex].pos + nodesDelta[duplicatedNodeNumber].GetNorm() * minEdgeLen * splitNodeBiasRatio);
        originalNodeIndices.push_back(nodeIndex);
      }
    }
    // all needed vertices are duplicated
    updatedGeomMesh->LoadGeom(mesh->vertices.data(), mesh->indices.data(),
      mesh->vertices.size(), mesh->GetCellsCount());
    updatedGeomMesh->BuildTopologyInfos();
    updatedGeomMesh->BuildUpdatedCellIndices(mesh->indices.data(), mesh->GetCellsCount());
  }
  void MarkGroupOfCells(IndexType cellIndex, std::map<IndexType, IndexType>& incidentCellsMarkers, IndexType cellsMarker);

  void GetIncidentCellsMarkers(IndexType nodeIndex, std::map<IndexType, IndexType>& incidentCellsMarkers, IndexType& cellMarkersCount)
  {
    usedCells.clear();
    // initialize dfs
    for (IndexType incidentCellNumber = 0; incidentCellNumber < geomMesh->nodeTopologyInfos[nodeIndex].incidentCellsCount; ++incidentCellNumber)
    {
      IndexType incidentCellIndex = geomMesh->nodeTopologyInfos[nodeIndex].incidentCells[incidentCellNumber];
      usedCells[incidentCellIndex] = false;
      incidentCellsMarkers[incidentCellIndex] = 0;
    }

    cellMarkersCount = 0;
    for (IndexType incidentCellNumber = 0; incidentCellNumber < geomMesh->nodeTopologyInfos[nodeIndex].incidentCellsCount; ++incidentCellNumber)
    {
      IndexType incidentCellIndex = geomMesh->nodeTopologyInfos[nodeIndex].incidentCells[incidentCellNumber];
      if (!usedCells[incidentCellIndex])
      {
        MarkGroupOfCells(incidentCellIndex, incidentCellsMarkers, cellMarkersCount);
        ++cellMarkersCount;
      }
    }
  }

  void BuildAdditionalTopology(GeomMesh<Space>* const geomMesh, MeshIO<Space>* const mesh);
};

template <>
void UnstructuredMeshBuilderCommon<Space2>::MarkGroupOfCells(IndexType cellIndex,
  std::map<IndexType, IndexType>& incidentCellsMarkers, IndexType cellsMarker)
{
  usedCells[cellIndex] = true;
  incidentCellsMarkers[cellIndex] = cellsMarker;

  for (IndexType faceNumber = 0; faceNumber < Space2::EdgesPerCell; ++faceNumber)
  {
    IndexType neighbourCellIndex = geomMesh->additionalCellInfos[cellIndex].neighbouringEdges[faceNumber].correspondingCellIndex;
    IndexType interactionType = geomMesh->additionalCellInfos[cellIndex].neighbouringEdges[faceNumber].interactionType;
    if (neighbourCellIndex != IndexType(-1) && !duplicateNodes && interactionType == 0)
    {
      std::map<IndexType, bool>::iterator it = usedCells.find(neighbourCellIndex);
      if (it != usedCells.end() && !it->second)
      {
        MarkGroupOfCells(neighbourCellIndex, incidentCellsMarkers, cellsMarker);
      }
    }
  }
}

template <>
void UnstructuredMeshBuilderCommon<Space3>::MarkGroupOfCells(IndexType cellIndex,
  std::map<IndexType, IndexType>& incidentCellsMarkers, IndexType cellsMarker)
{
  usedCells[cellIndex] = true;
  incidentCellsMarkers[cellIndex] = cellsMarker;

  for (IndexType faceNumber = 0; faceNumber < Space3::FacesPerCell; ++faceNumber)
  {
    IndexType neighbourCellIndex = geomMesh->additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
    IndexType interactionType    = geomMesh->additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;
    if (neighbourCellIndex != IndexType(-1) && !duplicateNodes && interactionType == 0)
    {
      std::map<IndexType, bool>::iterator it = usedCells.find(neighbourCellIndex);
      if (it != usedCells.end() && !it->second)
      {
        MarkGroupOfCells(neighbourCellIndex, incidentCellsMarkers, cellsMarker);
      }
    }
  }
}

template <>
void UnstructuredMeshBuilderCommon<Space2>::BuildAdditionalTopology(GeomMesh<Space2>* const geomMesh, 
  MeshIO<Space2>* const mesh)
{
  geomMesh->BuildAdditionalTopology(mesh->contactEdges.data(), mesh->contactEdgesCount.data(), mesh->contactTypesCount,
    mesh->boundaryEdges.data(), mesh->boundaryEdgesCount.data(), mesh->boundaryTypesCount);
}

template <>
void UnstructuredMeshBuilderCommon<Space3>::BuildAdditionalTopology(GeomMesh<Space3>* const geomMesh, 
  MeshIO<Space3>* const mesh)
{
  geomMesh->BuildAdditionalTopology(mesh->contactFaces.data(), mesh->contactFacesCount.data(), mesh->contactTypesCount, 
    mesh->boundaryFaces.data(), mesh->boundaryFacesCount.data(), mesh->boundaryTypesCount);
}

template <>
class UnstructuredMeshBuilder<Space2>: public UnstructuredMeshBuilderCommon<Space2>
{
public:
  SPACE2_TYPEDEFS
  typedef Space2 Space;

  UnstructuredMeshBuilder(bool duplicateNodes, Scalar splitNodeBiasRatio): 
    UnstructuredMeshBuilderCommon<Space>(duplicateNodes, splitNodeBiasRatio)
  {
  }

  virtual ~UnstructuredMeshBuilder() = default;

  typedef MeshIO<Space>::EdgePairIndices    EdgePairIndices;
  typedef MeshIO<Space>::BoundaryEdge       BoundaryEdge;
  typedef GeomMesh<Space>::EdgeLocation     EdgeLocation;
  typedef GeomMesh<Space>::EdgeIndices      EdgeIndices;
  typedef GeomMesh<Space>::EdgeLocationPair EdgeLocationPair;

protected:
  std::vector<BoundaryEdge*> boundaryData;
  std::vector<EdgePairIndices*> contactData;

  void BuildBoundary(MeshIO<Space2>* const mesh) override
  {
    std::vector<BoundaryEdge> boundaryEdges;
    IndexType offset = 0;
    IndexType boundariesCount = 0;

    for (IndexType phase = 0; phase < 2; ++phase)
    {
      if (phase == 1)
      {
        boundaryEdges.reserve(boundariesCount);
        offset = 0;
      }

      for (IndexType boundaryType = 0; boundaryType < mesh->boundaryTypesCount; ++boundaryType)
      {
        IndexType boundariesAdded = 0;
        for (IndexType boundaryNumber = 0; boundaryNumber < mesh->boundaryEdgesCount[boundaryType]; ++boundaryNumber)
        {
          IndexType edgeIndex = geomMesh->GetEdgeIndex(mesh->boundaryEdges[offset + boundaryNumber].nodeIndices);
          assert(edgeIndex != IndexType(-1));
          BoundaryEdge boundaryEdge = mesh->boundaryEdges[offset + boundaryNumber];
          std::sort(boundaryEdge.nodeIndices, boundaryEdge.nodeIndices + Space::NodesPerEdge);

          for (IndexType incidentCellNumber = 0; incidentCellNumber < geomMesh->edgeTopologyInfos[edgeIndex].incidentCellsCount; ++incidentCellNumber)
          {
            IndexType originalCellIndex = geomMesh->originalCellIndices[geomMesh->edgeTopologyInfos[edgeIndex].incidentCells[incidentCellNumber]];
            IndexType updatedCellIndex = updatedGeomMesh->updatedCellIndices[originalCellIndex];
            if (phase == 1)
            {
              BoundaryEdge edgeIndices;

              for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; ++edgeNumber)
              {
                IndexType edgeNodes[Space::NodesPerEdge];
                updatedGeomMesh->GetCellEdgeNodes(updatedCellIndex, edgeNumber, edgeNodes);
                IndexType originalEdgeNodes[Space::NodesPerEdge];
                for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
                {
                  originalEdgeNodes[nodeNumber] = originalNodeIndices[edgeNodes[nodeNumber]];
                }
                std::sort(originalEdgeNodes, originalEdgeNodes + Space::NodesPerEdge);
                bool found = true;
                for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
                {
                  if (originalEdgeNodes[nodeNumber] != boundaryEdge.nodeIndices[nodeNumber]) found = false;
                }
                if (found)
                {
                  std::copy(edgeNodes, edgeNodes + Space::NodesPerEdge, edgeIndices.nodeIndices);
                  break;
                }
              }
              boundaryEdges.push_back(edgeIndices);
            }
            boundariesAdded++;
          }
        }
        offset += mesh->boundaryEdgesCount[boundaryType];
        if (phase == 1)
        {
          mesh->boundaryEdgesCount[boundaryType] = boundariesAdded;
        }
        boundariesCount += boundariesAdded;
      }
    }
    mesh->boundaryEdges = boundaryEdges;
  }

  void BuildContacts(MeshIO<Space2>* const mesh) override
  {
    if (!duplicateNodes && mesh->contactTypesCount > 0)
    {
      mesh->contactEdges.erase(mesh->contactEdges.begin(), mesh->contactEdges.begin() + mesh->contactEdgesCount[0]);
      mesh->contactEdgesCount[0] = 0;
    }

    if (duplicateNodes)
    {
      if (mesh->contactTypesCount > 0 && mesh->contactEdgesCount[0] > 0)
      {
        std::cout << "Glue contacts have been already added" << std::endl;
      } else
      {
        typedef GeomMesh<Space>::EdgePairIndices EdgePairIndices;
        std::vector< EdgePairIndices > glueContacts;
        // add glue contacts
        for (IndexType cellIndex = 0; cellIndex < geomMesh->cells.size(); ++cellIndex)
        {
          for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; ++edgeNumber)
          {
            IndexType correspondingCellIndex = geomMesh->additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex;
            IndexType interactionType = geomMesh->additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;
            if (correspondingCellIndex != IndexType(-1) && interactionType == 0 && cellIndex < correspondingCellIndex)
            {
              EdgePairIndices edgePair;
              for (IndexType pair = 0; pair < 2; ++pair)
                geomMesh->GetCellEdgeNodes(cellIndex, edgeNumber, edgePair.edges[pair].nodeIndices);
              glueContacts.push_back(edgePair);
            }
          }
        }
        if (mesh->contactTypesCount == 0)
        {
          mesh->contactTypesCount++;
          mesh->contactEdgesCount.resize(1);
        }
        mesh->contactEdgesCount[0] = glueContacts.size();
        mesh->contactEdges.insert(mesh->contactEdges.begin(), glueContacts.begin(), glueContacts.end());
      }
    }

    IndexType offset = 0;
    for (IndexType contactType = 0; contactType < mesh->contactTypesCount; ++contactType)
    {
      for (IndexType contactNumber = 0; contactNumber < mesh->contactEdgesCount[contactType]; ++contactNumber)
      {
        if (duplicateNodes || contactType > 0)
        {
          EdgePairIndices contactEdge = mesh->contactEdges[offset + contactNumber];
          IndexType edgeIndex = geomMesh->GetEdgeIndex(contactEdge.edges[0].nodeIndices);
          assert(contactEdge.edges[0] == contactEdge.edges[1]);
          std::sort(contactEdge.edges[0].nodeIndices, contactEdge.edges[0].nodeIndices + Space::NodesPerEdge);

          EdgePairIndices splittedContactEdge;

          assert(edgeIndex != IndexType(-1));
          if (geomMesh->edgeTopologyInfos[edgeIndex].incidentCellsCount < 2) continue;
          for (IndexType incidentCellNumber = 0; incidentCellNumber < 2; ++incidentCellNumber)
          {
            IndexType originalCellIndex = geomMesh->originalCellIndices[geomMesh->edgeTopologyInfos[edgeIndex].incidentCells[incidentCellNumber]];
            IndexType updatedCellIndex = updatedGeomMesh->updatedCellIndices[originalCellIndex];
            for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; ++edgeNumber)
            {
              IndexType edgeNodes[Space::NodesPerEdge];
              updatedGeomMesh->GetCellEdgeNodes(updatedCellIndex, edgeNumber, edgeNodes);
              IndexType originalEdgeNodes[Space::NodesPerEdge];
              for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
              {
                originalEdgeNodes[nodeNumber] = originalNodeIndices[edgeNodes[nodeNumber]];
              }
              std::sort(originalEdgeNodes, originalEdgeNodes + Space::NodesPerEdge);
              bool found = true;
              for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; ++nodeNumber)
              {
                if (originalEdgeNodes[nodeNumber] != contactEdge.edges[0].nodeIndices[nodeNumber]) found = false;
              }
              if (found)
              {
                std::copy(edgeNodes, edgeNodes + Space::NodesPerEdge, splittedContactEdge.edges[incidentCellNumber].nodeIndices);
                break;
              }
            }
          }
          mesh->contactEdges[offset + contactNumber] = splittedContactEdge;
          assert(splittedContactEdge.edges[0] != splittedContactEdge.edges[1]);
        }
      }
      offset += mesh->contactEdgesCount[contactType];
    }
  }
};

template <>
class UnstructuredMeshBuilder<Space3>: public UnstructuredMeshBuilderCommon<Space3>
{
public:
  SPACE3_TYPEDEFS
  typedef Space3 Space;

  UnstructuredMeshBuilder(bool duplicateNodes, Scalar splitNodeBiasRatio): 
    UnstructuredMeshBuilderCommon<Space>(duplicateNodes, splitNodeBiasRatio)
  {
  }

  virtual ~UnstructuredMeshBuilder() = default;

  typedef MeshIO<Space>::FacePairIndices    FacePairIndices;
  typedef MeshIO<Space>::BoundaryFace       BoundaryFace;
  typedef GeomMesh<Space>::FaceLocation     FaceLocation;
  typedef GeomMesh<Space>::FaceIndices      FaceIndices;
  typedef GeomMesh<Space>::FaceLocationPair FaceLocationPair;

  using UnstructuredMeshBuilderCommon<Space>::duplicateNodes;
  using UnstructuredMeshBuilderCommon<Space>::splitNodeBiasRatio;

protected:
  void BuildBoundary(MeshIO<Space3>* const mesh) override
  {
    std::vector<BoundaryFace> boundaryFaces;
    IndexType offset = 0;
    IndexType boundariesCount = 0;

    for (IndexType phase = 0; phase < 2; ++phase)
    {
      if (phase == 1)
      {
        boundaryFaces.reserve(boundariesCount);
        offset = 0;
      }

      for (IndexType boundaryType = 0; boundaryType < mesh->boundaryTypesCount; ++boundaryType)
      {
        IndexType boundariesAdded = 0;
        for (IndexType boundaryNumber = 0; boundaryNumber < mesh->boundaryFacesCount[boundaryType]; ++boundaryNumber)
        {
          IndexType faceIndex = geomMesh->GetFaceIndex(mesh->boundaryFaces[offset + boundaryNumber].nodeIndices);
          assert(faceIndex != IndexType(-1));
          BoundaryFace boundaryFace = mesh->boundaryFaces[offset + boundaryNumber];
          std::sort(boundaryFace.nodeIndices, boundaryFace.nodeIndices + Space::NodesPerFace);

          for (IndexType incidentCellNumber = 0; incidentCellNumber < geomMesh->faceTopologyInfos[faceIndex].incidentCellsCount; ++incidentCellNumber)
          {
            IndexType originalCellIndex = geomMesh->originalCellIndices[geomMesh->faceTopologyInfos[faceIndex].incidentCells[incidentCellNumber]];
            IndexType updatedCellIndex = updatedGeomMesh->updatedCellIndices[originalCellIndex];
            if (phase == 1)
            {
              BoundaryFace faceIndices;

              for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
              {
                IndexType faceNodes[Space::NodesPerFace];
                updatedGeomMesh->GetCellFaceNodes(updatedCellIndex, faceNumber, faceNodes);
                IndexType originalFaceNodes[Space::NodesPerFace];
                for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
                {
                  originalFaceNodes[nodeNumber] = originalNodeIndices[faceNodes[nodeNumber]];
                }
                std::sort(originalFaceNodes, originalFaceNodes + Space::NodesPerFace);
                bool found = true;
                for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
                {
                  if (originalFaceNodes[nodeNumber] != boundaryFace.nodeIndices[nodeNumber]) found = false;
                }
                if (found)
                {
                  std::copy(faceNodes, faceNodes + Space::NodesPerFace, faceIndices.nodeIndices);
                  break;
                }
              }
              boundaryFaces.push_back(faceIndices);
            }
            boundariesAdded++; 
          }
        }
        offset += mesh->boundaryFacesCount[boundaryType];
        if (phase == 1)
        {
          mesh->boundaryFacesCount[boundaryType] = boundariesAdded;
        }
        boundariesCount += boundariesAdded;
      }
    }
    mesh->boundaryFaces = boundaryFaces;
  }

  bool RegulateFaceNodeIndices(IndexType faceNodes[Space::NodesPerFace], 
    IndexType testOriginalNodes[Space::NodesPerFace], IndexType testUpdatedNodes[Space::NodesPerFace],
    IndexType* regulatedNodes) const
  {
    for (IndexType i = 0; i + 1 < Space::NodesPerFace; ++i)
    {
      for (IndexType j = i + 1; j < Space::NodesPerFace; ++j)
      {
        if (testOriginalNodes[i] > testOriginalNodes[j])
        {
          std::swap(testOriginalNodes[i], testOriginalNodes[j]);
          std::swap(testUpdatedNodes[i], testUpdatedNodes[j]);
        }
      }
    }
    bool found = true;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
    {
      if (faceNodes[nodeNumber] != testOriginalNodes[nodeNumber]) found = false;
      regulatedNodes[nodeNumber] = testUpdatedNodes[nodeNumber];
    }
    return found;
  }

  void BuildContacts(MeshIO<Space>* const mesh) override
  {
    if (!duplicateNodes && mesh->contactTypesCount > 0)
    {
      mesh->contactFaces.erase(mesh->contactFaces.begin(), mesh->contactFaces.begin() + mesh->contactFacesCount[0]);
      mesh->contactFacesCount[0] = 0;
    } 

    if (duplicateNodes)
    {
      if (mesh->contactTypesCount > 0 && mesh->contactFacesCount[0] > 0)
      {
        std::cout << "Glue contacts have been already added" << std::endl;
      } else
      {
        typedef GeomMesh<Space>::FacePairIndices FacePairIndices;
        std::vector< FacePairIndices > glueContacts;
        // add glue contacts
        for (IndexType cellIndex = 0; cellIndex < geomMesh->cells.size(); ++cellIndex)
        {
          for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
          {
            IndexType correspondingCellIndex = geomMesh->additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
            IndexType interactionType = geomMesh->additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;
            if (correspondingCellIndex != IndexType(-1) &&
                interactionType == 0 && cellIndex < correspondingCellIndex)
            {
              FacePairIndices facePair;
              for (IndexType pair = 0; pair < 2; ++pair)
                geomMesh->GetCellFaceNodes(cellIndex, faceNumber, facePair.faces[pair].nodeIndices);
              glueContacts.push_back(facePair);
            }
          }
        }
        if (mesh->contactTypesCount == 0)
        {
          mesh->contactTypesCount++;
          mesh->contactFacesCount.resize(1);
        }
        mesh->contactFacesCount[0] = glueContacts.size();
        mesh->contactFaces.insert(mesh->contactFaces.begin(), glueContacts.begin(), glueContacts.end());
      }
    }

    IndexType offset = 0;
    for (IndexType contactType = 0; contactType < mesh->contactTypesCount; ++contactType)
    {
      for (IndexType contactNumber = 0; contactNumber < mesh->contactFacesCount[contactType]; ++contactNumber)
      {
        if (duplicateNodes || contactType > 0)
        {
          FacePairIndices contactFace = mesh->contactFaces[offset + contactNumber];
          IndexType faceIndex = geomMesh->GetFaceIndex(contactFace.faces[0].nodeIndices);
          assert(contactFace.faces[0] == contactFace.faces[1]);
          std::sort(contactFace.faces[0].nodeIndices, contactFace.faces[0].nodeIndices + Space::NodesPerFace);

          FacePairIndices splittedContactFace;

          assert(faceIndex != IndexType(-1));
          if (geomMesh->faceTopologyInfos[faceIndex].incidentCellsCount < 2) continue;
          for (IndexType incidentCellNumber = 0; incidentCellNumber < 2; ++incidentCellNumber)
          {
            IndexType originalCellIndex = geomMesh->originalCellIndices[geomMesh->faceTopologyInfos[faceIndex].incidentCells[incidentCellNumber]];
            IndexType updatedCellIndex = updatedGeomMesh->updatedCellIndices[originalCellIndex];
            bool found = false;
            for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
            {
              IndexType faceNodes[Space::NodesPerFace];
              updatedGeomMesh->GetCellFaceNodes(updatedCellIndex, faceNumber, faceNodes);
              IndexType originalFaceNodes[Space::NodesPerFace];
              for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
              {
                originalFaceNodes[nodeNumber] = originalNodeIndices[faceNodes[nodeNumber]];
              }
              if (RegulateFaceNodeIndices(contactFace.faces[0].nodeIndices, originalFaceNodes, faceNodes,
                splittedContactFace.faces[incidentCellNumber].nodeIndices))
              {
                found = true;
                break;
              }
            }
            assert(found);
          }
          mesh->contactFaces[offset + contactNumber] = splittedContactFace;
          assert(splittedContactFace.faces[0] != splittedContactFace.faces[1]);
        }
      }
      offset += mesh->contactFacesCount[contactType];
    }
  }
};

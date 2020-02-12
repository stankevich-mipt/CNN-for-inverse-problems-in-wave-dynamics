#pragma once

#include "../GeomMesh/MeshIO/Distributed/DistributedMeshIO.h"
#include "../GeomMesh/GeomMesh/GeomMesh.h"

#include <algorithm>

template <typename Space>
class MeshChecker
{
public:
  SPACE_TYPEDEFS

  MeshChecker(const std::string& meshBaseName, IndexType domainsCount):
    minEdgeLen(std::numeric_limits<Scalar>::max()), domains(domainsCount, DistributedMeshIO<Space>(domainsCount)), meshBaseName(meshBaseName)
  {
    std::cout << "Checking mesh..." << std::endl;
  }

  void Check()
  {
    LoadDomains();
    // contacts
    std::cout << "  Whole mesh... " << (CheckContacts(baseMesh) ? "Ok" : "Error") << std::endl;
    for (IndexType domainIndex = 0; domainIndex < domains.size(); ++domainIndex)
    {
      std::cout << "  Domain " << domainIndex << "... " << (CheckContacts(domains[domainIndex]) ? "Ok" : "Error") << std::endl;
    }
    
    std::cout << "  Transition info... " << (CheckTransitionInfo() ? "Ok" : "Error") << std::endl;
  }

public:
  void LoadDomains()
  {
    for (IndexType domainIndex = 0; domainIndex < domains.size(); ++domainIndex)
    {
      char domainString[256];
      sprintf(domainString, "%.3lu", domainIndex);

      std::string meshName = meshBaseName;
      ReplaceSubstring(meshName, "<domain>", std::string(domainString));
      domains[domainIndex].Load(AddExtensionToFileName(meshName, ".mesh"));
    }
    std::string meshName = meshBaseName;
    ReplaceSubstring(meshName, "<domain>", "");
    baseMesh.Load(AddExtensionToFileName(meshName, ".mesh"));
    
    for (IndexType cellIndex = 0; cellIndex < baseMesh.GetCellsCount(); ++cellIndex)
    {
      Vector points[Space::NodesPerCell];
      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
      {
        points[nodeNumber] = baseMesh.vertices[baseMesh.indices[Space::NodesPerCell * cellIndex + nodeNumber]];
      }
      minEdgeLen = Min(minEdgeLen, GetMinEdgeLen(points));
    }
  }

  Scalar GetMinEdgeLen(Vector points[Space::NodesPerCell]) const
  {
    Scalar minDist = std::numeric_limits<Scalar>::max();
    for (IndexType firstPointNumber = 0; firstPointNumber + 1 < Space::NodesPerCell; ++firstPointNumber)
    {
      for (IndexType secondPointNumber = firstPointNumber + 1; secondPointNumber < Space::NodesPerCell; ++secondPointNumber)
      {
        minDist = Min(minDist, (points[secondPointNumber] - points[firstPointNumber]).Len());
      }
    }
    return minDist;
  }

  bool CheckContacts(const MeshIO<Space>& mesh);

  bool CheckTransitionInfo()
  {
    bool isOk = true;
    std::vector< GeomMesh<Space> > geomMeshs(domains.size());
    for (IndexType domainIndex = 0; domainIndex < domains.size(); ++domainIndex)
    {
      geomMeshs[domainIndex].LoadGeom(domains[domainIndex].vertices.data(), domains[domainIndex].indices.data(),
        domains[domainIndex].vertices.size(), domains[domainIndex].GetCellsCount());
      geomMeshs[domainIndex].BuildTopologyInfos();
    }

    for (IndexType sendingDomainIndex = 0; sendingDomainIndex < domains.size(); ++sendingDomainIndex)
    {
      for (IndexType receivingDomainIndex = 0; receivingDomainIndex < domains.size(); ++receivingDomainIndex)
      {
        TransitionInfo<Space>& sendingInfo = domains[sendingDomainIndex].transitionInfos[receivingDomainIndex];        
        // check node indices
        std::map<IndexType, IndexType> transitionNodeMap;
        for (IndexType nodeIndex = 0; nodeIndex < sendingInfo.transitionNodes.size(); ++nodeIndex)
        {
          transitionNodeMap[sendingInfo.transitionNodes[nodeIndex].nativeIndex] = sendingInfo.transitionNodes[nodeIndex].targetIndex;
        }

        for (IndexType nodeIndex = 0; nodeIndex < sendingInfo.nodesIndices.size(); ++nodeIndex)
        {
          Vector sNode = domains[sendingDomainIndex].vertices[sendingInfo.nodesIndices[nodeIndex]];
          Vector rNode = domains[receivingDomainIndex].vertices[transitionNodeMap[sendingInfo.nodesIndices[nodeIndex]]];
          if ((sNode - rNode).Len() > Scalar(1e-1) * minEdgeLen)
          {
            isOk = false;
          }
        }

        // check cell indices
        for (IndexType cellIndex = 0; cellIndex < sendingInfo.cells.size(); ++cellIndex)
        {
          IndexType translatedIndices[Space::NodesPerCell];
          for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
          {
            translatedIndices[nodeNumber] = transitionNodeMap[sendingInfo.cells[cellIndex].incidentNodes[nodeNumber]];
          }
          bool ordered = true;
          for (IndexType nodeNumber = 0; nodeNumber + 1 < Space::NodesPerCell; ++nodeNumber)
          {
            if (sendingInfo.cells[cellIndex].incidentNodes[nodeNumber] > sendingInfo.cells[cellIndex].incidentNodes[nodeNumber + 1]) ordered = false;
            if (translatedIndices[nodeNumber] > translatedIndices[nodeNumber + 1]) ordered = false;
          }
          if (!ordered)
          {
            isOk = false;
            std::cout << "Bad cell orientation\n";
          }         
        }
      }
    }
    return isOk;
  }

  Scalar minEdgeLen;
  std::vector< DistributedMeshIO<Space> > domains;
  MeshIO<Space> baseMesh;
  std::string meshBaseName;
};

template <>
bool MeshChecker<Space2>::CheckContacts(const MeshIO<Space2>& mesh)
{
  typedef GeomMesh<Space2>::EdgePairIndices EdgePairIndices;
  bool isOk = true;

  IndexType offset = 0;
  for (IndexType contactTypesIndex = 0; contactTypesIndex < mesh.contactTypesCount; ++contactTypesIndex)
  {
    for (IndexType contactEdgeIndex = 0; contactEdgeIndex < mesh.contactEdgesCount[contactTypesIndex]; ++contactEdgeIndex)
    {
      EdgePairIndices edgePairIndices = mesh.contactEdges[offset + contactEdgeIndex];
      if (edgePairIndices.edges[0] != edgePairIndices.edges[1])
      {
        Vector leftSidePoints[2];
        Vector rightSidePoints[2];
        for (IndexType nodeNumber = 0; nodeNumber < 2; ++nodeNumber)
        {
          leftSidePoints[nodeNumber] = mesh.vertices[edgePairIndices.edges[0].nodeIndices[nodeNumber]];
          rightSidePoints[nodeNumber] = mesh.vertices[edgePairIndices.edges[1].nodeIndices[nodeNumber]];
        }

        Scalar eps = minEdgeLen * Scalar(1e-1);

        if ((leftSidePoints[0] - rightSidePoints[0]).Len() < eps && (leftSidePoints[1] - rightSidePoints[1]).Len() < eps) continue;
        if ((leftSidePoints[0] - rightSidePoints[1]).Len() < eps && (leftSidePoints[1] - rightSidePoints[0]).Len() < eps) continue;

        std::cout << "Nonconform contact: "
          << "contactTypeIndex " << contactTypesIndex << ";"
          << "contactEdgeIndex " << contactEdgeIndex << ";";
        std::cout << "indices: ";
        for (IndexType edgeNumber = 0; edgeNumber < 2; ++edgeNumber)
        {
          std::cout << edgePairIndices.edges[edgeNumber] << ";";
        }
        std::cout << std::endl;
        isOk = false;
      } else
      {
        std::cout << "Nonduplicated contact: " 
          << "contactTypeIndex " << contactTypesIndex << ";"
          << "contactEdgeIndex " << contactEdgeIndex << ";";

        for (IndexType edgeNumber = 0; edgeNumber < 2; ++edgeNumber)
        {
          std::cout << edgePairIndices.edges[edgeNumber] << ";";
        }
        std::cout << std::endl;
        isOk = false;
      }
    }
    offset += mesh.contactEdgesCount[contactTypesIndex];
  }

  return isOk;
}

template <>
bool MeshChecker<Space3>::CheckContacts(const MeshIO<Space3>& mesh)
{
  typedef GeomMesh<Space3>::FacePairIndices FacePairIndices;
  bool isOk = true;

  IndexType offset = 0;
  for (IndexType contactTypesIndex = 0; contactTypesIndex < mesh.contactTypesCount; ++contactTypesIndex)
  {
    for (IndexType contactFaceIndex = 0; contactFaceIndex < mesh.contactFacesCount[contactTypesIndex]; ++contactFaceIndex)
    {
      FacePairIndices facePairIndices = mesh.contactFaces[offset + contactFaceIndex];
      if (facePairIndices.faces[0] != facePairIndices.faces[1])
      {
        Vector leftSidePoints[Space3::NodesPerFace];
        Vector rightSidePoints[Space3::NodesPerFace];
        for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
        {
          leftSidePoints[nodeNumber] = mesh.vertices[facePairIndices.faces[0].nodeIndices[nodeNumber]];
          rightSidePoints[nodeNumber] = mesh.vertices[facePairIndices.faces[1].nodeIndices[nodeNumber]];
        }

        Scalar eps = 1.0 / 5.0;

        if ((leftSidePoints[0] - rightSidePoints[0]).Len() < eps && 
            (leftSidePoints[1] - rightSidePoints[1]).Len() < eps &&
            (leftSidePoints[2] - rightSidePoints[2]).Len() < eps) continue;

        std::cout << "Nonconform contact: "
          << "contactTypeIndex " << contactTypesIndex << ";"
          << "contactEdgeIndex " << contactFaceIndex << ";";
        std::cout << "indices: ";
        for (IndexType faceNumber = 0; faceNumber < 2; ++faceNumber)
        {
          std::cout << facePairIndices.faces[faceNumber] << ";";
        }
        std::cout << std::endl;
        isOk = false;
      }
      else
      {
        std::cout << "Nonduplicated contact: "
          << "contactTypeIndex " << contactTypesIndex << ";"
          << "contactEdgeIndex " << contactFaceIndex << ";";

        for (IndexType faceNumber = 0; faceNumber < 2; ++faceNumber)
        {
          std::cout << facePairIndices.faces[faceNumber] << ";";
        }
        std::cout << std::endl;
        isOk = false;
      }
    }
    offset += mesh.contactFacesCount[contactTypesIndex];
  }

  return isOk;
}

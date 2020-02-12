#include <algorithm>
#include <vector>
#include "../../Utils/Utils.h"

template<typename T, typename IndexType>
class DuplicateRemover
{
public:
  void AddElement(const T& elem)
  {
    elements.push_back(elem);
  }
  IndexType RemoveDuplicates()
  {
    return ::RemoveDuplicates(&elements);
  }
  T& operator [](IndexType num)
  {
    return elements[num];
  }
  T* GetElements()
  {
    return elements.data();
  }
  IndexType GetElementsCount()
  {
    return IndexType(elements.size());
  }
private:
  std::vector<T> elements;
};

template<typename IndexType>
struct TopologyReconstructor
{
private:
  struct NodeNeighbours
  {
    IndexType incidentCellsCount;
    IndexType incidentEdgesCount;
    IndexType incidentFacesCount;

    IndexType incidentCellsStart;
    IndexType incidentEdgesStart;
    IndexType incidentFacesStart;
  };

  struct EdgeNeighbours
  {
    EdgeNeighbours()
    {
      incidentCellsCount = 0;
    }
    IndexType incidentCellsCount;
    IndexType incidentCellsStart;
  };

  struct EdgeIndices
  {
    EdgeIndices()
    {
    }
    EdgeIndices(IndexType node0, IndexType node1)
    {
      nodes[0] = node0;
      nodes[1] = node1;
      std::sort(nodes, nodes + 2);
    }
    bool operator ==(const EdgeIndices &e) const
    {
      return ((nodes[0] == e.nodes[0]) && (nodes[1] == e.nodes[1]));
    }
    bool operator < (const EdgeIndices &e) const
    {
      if(nodes[0] < e.nodes[0]) return 1;
      if((nodes[0] == e.nodes[0]) && (nodes[1] < e.nodes[1])) return 1;
      return 0;
    }
    IndexType nodes[2];
  };

  struct FaceNeighbours
  {
    FaceNeighbours()
    {
      incidentCellsCount = 0;
    }

    IndexType incidentCellsCount;
    IndexType incidentCellsStart;
  };

  struct FaceIndices
  {
    FaceIndices()
    {
    }
    FaceIndices(IndexType node0, IndexType node1, IndexType node2)
    {
      nodes[0] = node0;
      nodes[1] = node1;
      nodes[2] = node2;
      std::sort(nodes, nodes + 3);
    }

    bool operator ==(const FaceIndices &f) const
    {
      return ((nodes[0] == f.nodes[0]) && (nodes[1] == f.nodes[1]) && (nodes[2] == f.nodes[2]));
    }
    bool operator < (const FaceIndices &f) const
    {
      if(nodes[0] < f.nodes[0]) return 1;
      if((nodes[0] == f.nodes[0]) && (nodes[1] < f.nodes[1])) return 1;
      if((nodes[0] == f.nodes[0]) && (nodes[1] == f.nodes[1]) && (nodes[2] < f.nodes[2])) return 1;
      return 0;
    }

    IndexType nodes[3];
  };

  struct CellIndices
  {
    CellIndices()
    {
    }
    CellIndices(IndexType node0, IndexType node1, IndexType node2, IndexType node3)
    {
      nodes[0] = node0;
      nodes[1] = node1;
      nodes[2] = node2;
      nodes[3] = node3;
      std::sort(nodes, nodes + 4);
    }
    bool operator < (const CellIndices &c) const
    {
      if(nodes[0] < c.nodes[0]) return 1;
      if((nodes[0] == c.nodes[0]) && (nodes[1] < c.nodes[1])) return 1;
      if((nodes[0] == c.nodes[0]) && (nodes[1] == c.nodes[1]) && (nodes[2] < c.nodes[2])) return 1;
      if((nodes[0] == c.nodes[0]) && (nodes[1] == c.nodes[1]) && (nodes[2] == c.nodes[2]) && (nodes[3] < c.nodes[3])) return 1;
      return 0;
    }
    IndexType nodes[4];
  };

  struct CellNeighbours
  {
    CellNeighbours()
    {
    }
    IndexType edges[6];
    IndexType faces[4];
  };
  void BuildCellIndices(IndexType *indexData, IndexType indices)
  {
    //building CellIndices and calculating total nodes count
    cellsCount = indices;
    cellIndices.reserve(indices);

    IndexType maxNodeIndex;

    if(indices > 0)
    {
      maxNodeIndex = indexData[0];
    }else
    {
      nodesCount = 0;
      return;
    }   
//    int s = cellIndices.size();

    for(IndexType k = 0; k < indices; k++)
    {
      cellIndices.push_back(CellIndices(indexData[4 * k + 0], indexData[4 * k + 1], indexData[4 * k + 2], indexData[4 * k + 3]));
      for(IndexType i = 0; i < 4; i++)
      {
        if(cellIndices[k].nodes[i] > maxNodeIndex)
        {
          maxNodeIndex = cellIndices[k].nodes[i];
        }
      }
    }
    std::sort(cellIndices.begin(), cellIndices.end());
    nodesCount = maxNodeIndex + 1; //number of nodes is more than max node index by one
  }
  void BuildEdgeIndices()
  {
    //building edge indices. add one by one each tetrahedra's edge and then remove duplicates
    for(IndexType i = 0; i < cellsCount; i++)
    {
      edgeIndices.AddElement(EdgeIndices(cellIndices[i].nodes[0], cellIndices[i].nodes[1]));
      edgeIndices.AddElement(EdgeIndices(cellIndices[i].nodes[1], cellIndices[i].nodes[2]));
      edgeIndices.AddElement(EdgeIndices(cellIndices[i].nodes[2], cellIndices[i].nodes[0]));
      edgeIndices.AddElement(EdgeIndices(cellIndices[i].nodes[0], cellIndices[i].nodes[3]));
      edgeIndices.AddElement(EdgeIndices(cellIndices[i].nodes[1], cellIndices[i].nodes[3]));
      edgeIndices.AddElement(EdgeIndices(cellIndices[i].nodes[2], cellIndices[i].nodes[3]));
    }
    edgesCount = edgeIndices.RemoveDuplicates();
  }
  void BuildFaceIndices()
  {
    //building face indices. just like edges
    for(IndexType i = 0; i < cellsCount; i++)
    {
      faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[0], cellIndices[i].nodes[1], cellIndices[i].nodes[2]));
      faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[0], cellIndices[i].nodes[1], cellIndices[i].nodes[3]));
      faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[0], cellIndices[i].nodes[2], cellIndices[i].nodes[3]));
      faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[1], cellIndices[i].nodes[2], cellIndices[i].nodes[3]));
    }
    facesCount = faceIndices.RemoveDuplicates();
  }
  void InitNeighbouringStructures()
  {
    //initialize topology structures
    nodeNeighbours.resize(nodesCount);
    for(IndexType i = 0; i < nodesCount; i++)
    {
      nodeNeighbours[i].incidentCellsCount = 0;
      nodeNeighbours[i].incidentEdgesCount = 0;
      nodeNeighbours[i].incidentFacesCount = 0;
    }

    edgeNeighbours.resize(edgesCount);
    for(IndexType i = 0; i < edgesCount; i++)
    {
      edgeNeighbours[i].incidentCellsCount = 0;
    }

    faceNeighbours.resize(facesCount);
    for(IndexType i = 0; i < facesCount; i++)
    {
      faceNeighbours[i].incidentCellsCount = 0;
    }

    cellNeighbours.resize(cellsCount);
  }
  void BuildNodeIncidentCells()
  {
    //calculate each node's incident(neighbouring) cells count and form global pool of cells incident to nodes
    totalNodeIncidentCells = 0;
    for(IndexType i = 0; i < cellsCount; i++)
    {
      nodeNeighbours[cellIndices[i].nodes[0]].incidentCellsCount++;
      nodeNeighbours[cellIndices[i].nodes[1]].incidentCellsCount++;
      nodeNeighbours[cellIndices[i].nodes[2]].incidentCellsCount++;
      nodeNeighbours[cellIndices[i].nodes[3]].incidentCellsCount++;
      //each cell increases the total number of incident cells by 4, because it is incident to 4 nodes.
      totalNodeIncidentCells += 4;
    }
    nodeIncidentCellsPool.resize(totalNodeIncidentCells);

    //assign a place in pool for each node to store its incident cells
    IndexType poolState = 0;
    for(IndexType i = 0; i < nodesCount; i++)
    {
      nodeNeighbours[i].incidentCellsStart = poolState;
      poolState += nodeNeighbours[i].incidentCellsCount;
      nodeNeighbours[i].incidentCellsCount = 0;
    }
    if(poolState != totalNodeIncidentCells)
    {
      //fail
    }


    //assign cell indices to global pool
    for(IndexType i = 0; i < cellsCount; i++)
    {
      IndexType node0 = cellIndices[i].nodes[0];
      nodeIncidentCellsPool[nodeNeighbours[node0].incidentCellsStart + nodeNeighbours[node0].incidentCellsCount++] = i;
      IndexType node1 = cellIndices[i].nodes[1];
      nodeIncidentCellsPool[nodeNeighbours[node1].incidentCellsStart + nodeNeighbours[node1].incidentCellsCount++] = i;
      IndexType node2 = cellIndices[i].nodes[2];
      nodeIncidentCellsPool[nodeNeighbours[node2].incidentCellsStart + nodeNeighbours[node2].incidentCellsCount++] = i;
      IndexType node3 = cellIndices[i].nodes[3];
      nodeIncidentCellsPool[nodeNeighbours[node3].incidentCellsStart + nodeNeighbours[node3].incidentCellsCount++] = i;
    }
  }
  void BuildNodeIncidentFaces()
  {
    //calculate each node's incident faces count and form global pool of faces incident to nodes. just as cells
    totalNodeIncidentFaces = 0;
    for(IndexType i = 0; i < facesCount; i++)
    {
      nodeNeighbours[faceIndices[i].nodes[0]].incidentFacesCount++;
      nodeNeighbours[faceIndices[i].nodes[1]].incidentFacesCount++;
      nodeNeighbours[faceIndices[i].nodes[2]].incidentFacesCount++;
      //each edge increases the total number of incident edges by 2, because it is incident to 2 nodes.
      totalNodeIncidentFaces += 3;
    }
    nodeIncidentFacesPool.resize(totalNodeIncidentFaces);

    //assign a place in pool for each node to store its incident cells
    IndexType poolState = 0;
    for(IndexType i = 0; i < nodesCount; i++)
    {
      nodeNeighbours[i].incidentFacesStart = poolState;
      poolState += nodeNeighbours[i].incidentFacesCount;
      nodeNeighbours[i].incidentFacesCount = 0;
    }
    if(poolState != totalNodeIncidentFaces)
    {
      //fail
    }


    //assign face indices to global pool
    for(IndexType i = 0; i < facesCount; i++)
    {
      IndexType node0 = faceIndices[i].nodes[0];
      nodeIncidentFacesPool[nodeNeighbours[node0].incidentFacesStart + nodeNeighbours[node0].incidentFacesCount++] = i;
      IndexType node1 = faceIndices[i].nodes[1];
      nodeIncidentFacesPool[nodeNeighbours[node1].incidentFacesStart + nodeNeighbours[node1].incidentFacesCount++] = i;
      IndexType node2 = faceIndices[i].nodes[2];
      nodeIncidentFacesPool[nodeNeighbours[node2].incidentFacesStart + nodeNeighbours[node2].incidentFacesCount++] = i;
    }
  }
  void BuildNodeIncidentEdges()
  {
    //calculate each node's incident edges count and form global pool of edges incident to nodes. just as faces and cells
    totalNodeIncidentEdges = 0;
    for(IndexType i = 0; i < edgesCount; i++)
    {
      nodeNeighbours[edgeIndices[i].nodes[0]].incidentEdgesCount++;
      nodeNeighbours[edgeIndices[i].nodes[1]].incidentEdgesCount++;
      //each edge increases the total number of incident edges by 2, because it is incident to 2 nodes.
      totalNodeIncidentEdges += 2;
    }

    nodeIncidentEdgesPool.resize(totalNodeIncidentEdges);

    //assign a place in pool for each node to store its incident cells
    IndexType poolState = 0;
    for(IndexType i = 0; i < nodesCount; i++)
    {
      nodeNeighbours[i].incidentEdgesStart = poolState;
      poolState += nodeNeighbours[i].incidentEdgesCount;
      nodeNeighbours[i].incidentEdgesCount = 0;
    }
    if(poolState != totalNodeIncidentEdges)
    {
      //fail
    }


    //assign edge indices to global pool
    for(IndexType i = 0; i < edgesCount; i++)
    {
      IndexType node0 = edgeIndices[i].nodes[0];
      nodeIncidentEdgesPool[nodeNeighbours[node0].incidentEdgesStart + nodeNeighbours[node0].incidentEdgesCount++] = i;
      IndexType node1 = edgeIndices[i].nodes[1];
      nodeIncidentEdgesPool[nodeNeighbours[node1].incidentEdgesStart + nodeNeighbours[node1].incidentEdgesCount++] = i;
    }
  }

  void BuildCellIncidentEdges_old()
  {
    for(IndexType i = 0; i < cellsCount; i++)
    {
      IndexType node0 = cellIndices[i].nodes[0];
      IndexType node1 = cellIndices[i].nodes[1];
      IndexType node2 = cellIndices[i].nodes[2];
      IndexType node3 = cellIndices[i].nodes[3];

      IndexType edge0 = -1;
      for(IndexType j = 0; j < nodeNeighbours[node0].incidentEdgesCount; j++)
      {
        IndexType edgeIndex = nodeIncidentEdgesPool[nodeNeighbours[node0].incidentEdgesStart + j];
        if((edgeIndices[edgeIndex].nodes[0] == node0) && (edgeIndices[edgeIndex].nodes[1] == node1))
        {
          edge0 = edgeIndex;
          break;
        }
      }
      if(edge0 == IndexType(-1))
      {
        //fail
      }else
        cellNeighbours[i].edges[0] = edge0;

      IndexType edge1 = IndexType(-1);
      for(IndexType j = 0; j < nodeNeighbours[node0].incidentEdgesCount; j++)
      {
        IndexType edgeIndex = nodeIncidentEdgesPool[nodeNeighbours[node0].incidentEdgesStart + j];
        if((edgeIndices[edgeIndex].nodes[0] == node0) && (edgeIndices[edgeIndex].nodes[1] == node2))
        {
          edge1 = edgeIndex;
          break;
        }
      }
      if(edge1 == IndexType(-1))
      {
        //fail
      }else
        cellNeighbours[i].edges[1] = edge1;

      IndexType edge2 = IndexType(-1);
      for(IndexType j = 0; j < nodeNeighbours[node0].incidentEdgesCount; j++)
      {
        IndexType edgeIndex = nodeIncidentEdgesPool[nodeNeighbours[node0].incidentEdgesStart + j];
        if((edgeIndices[edgeIndex].nodes[0] == node0) && (edgeIndices[edgeIndex].nodes[1] == node3))
        {
          edge2 = edgeIndex;
          break;
        }
      }
      if(edge2 == IndexType(-1))
      {
        //fail
      }else
        cellNeighbours[i].edges[2] = edge2;


      IndexType edge3 = IndexType(-1);
      for(IndexType j = 0; j < nodeNeighbours[node1].incidentEdgesCount; j++)
      {
        IndexType edgeIndex = nodeIncidentEdgesPool[nodeNeighbours[node1].incidentEdgesStart + j];
        if((edgeIndices[edgeIndex].nodes[0] == node1) && (edgeIndices[edgeIndex].nodes[1] == node2))
        {
          edge3 = edgeIndex;
          break;
        }
      }
      if(edge3 == IndexType(-1))
      {
        //fail
      }else
        cellNeighbours[i].edges[3] = edge3;

      IndexType edge4 = IndexType(-1);
      for(IndexType j = 0; j < nodeNeighbours[node1].incidentEdgesCount; j++)
      {
        IndexType edgeIndex = nodeIncidentEdgesPool[nodeNeighbours[node1].incidentEdgesStart + j];
        if((edgeIndices[edgeIndex].nodes[0] == node1) && (edgeIndices[edgeIndex].nodes[1] == node3))
        {
          edge4 = edgeIndex;
          break;
        }
      }
      if(edge4 == IndexType(-1))
      {
        //fail
      }else
        cellNeighbours[i].edges[4] = edge4;

      IndexType edge5 = IndexType(-1);
      for(IndexType j = 0; j < nodeNeighbours[node2].incidentEdgesCount; j++)
      {
        IndexType edgeIndex = nodeIncidentEdgesPool[nodeNeighbours[node2].incidentEdgesStart + j];
        if((edgeIndices[edgeIndex].nodes[0] == node2) && (edgeIndices[edgeIndex].nodes[1] == node3))
        {
          edge5 = edgeIndex;
          break;
        }
      }
      if(edge5 == IndexType(-1))
      {
        //fail
      }else
        cellNeighbours[i].edges[5] = edge5;
    }
  }


  void BuildCellIncidentEdges()
  {
    for(IndexType i = 0; i < edgesCount; i++)
    {
      IndexType edgeNode0 = edgeIndices[i].nodes[0];
      IndexType edgeNode1 = edgeIndices[i].nodes[1];

      for(IndexType j = 0; j < nodeNeighbours[edgeNode0].incidentCellsCount; j++)
      {
        IndexType cellIndex = nodeIncidentCellsPool[nodeNeighbours[edgeNode0].incidentCellsStart + j];

        IndexType cellNode0 = cellIndices[cellIndex].nodes[0];
        IndexType cellNode1 = cellIndices[cellIndex].nodes[1];
        IndexType cellNode2 = cellIndices[cellIndex].nodes[2];
        IndexType cellNode3 = cellIndices[cellIndex].nodes[3];

        IndexType edgeInCellIndex = -1;

        if((cellNode0 == edgeNode0) && (cellNode1 == edgeNode1))
          edgeInCellIndex = 0;
        if((cellNode0 == edgeNode0) && (cellNode2 == edgeNode1))
          edgeInCellIndex = 1;
        if((cellNode0 == edgeNode0) && (cellNode3 == edgeNode1))
          edgeInCellIndex = 2;
        if((cellNode1 == edgeNode0) && (cellNode2 == edgeNode1))
          edgeInCellIndex = 3;
        if((cellNode1 == edgeNode0) && (cellNode3 == edgeNode1))
          edgeInCellIndex = 4;
        if((cellNode2 == edgeNode0) && (cellNode3 == edgeNode1))
          edgeInCellIndex = 5;


        if(edgeInCellIndex != IndexType(-1))
          cellNeighbours[cellIndex].edges[edgeInCellIndex] = i;
      }
    }
  }
  void BuildCellIncidentFaces_old()
  {
    for(IndexType i = 0; i < cellsCount; i++)
    {
      IndexType node0 = cellIndices[i].nodes[0];
      IndexType node1 = cellIndices[i].nodes[1];
      IndexType node2 = cellIndices[i].nodes[2];
      IndexType node3 = cellIndices[i].nodes[3];

      for(IndexType j = 0; j < nodeNeighbours[node0].incidentFacesCount; j++)
      {
        IndexType faceIndex = nodeIncidentFacesPool[nodeNeighbours[node0].incidentFacesStart + j];

        faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[0], cellIndices[i].nodes[1], cellIndices[i].nodes[2]));
        faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[0], cellIndices[i].nodes[1], cellIndices[i].nodes[3]));
        faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[1], cellIndices[i].nodes[2], cellIndices[i].nodes[3]));
        faceIndices.AddElement(FaceIndices(cellIndices[i].nodes[2], cellIndices[i].nodes[0], cellIndices[i].nodes[3]));

        if((faceIndices[faceIndex].nodes[0] == node0) &&
           (faceIndices[faceIndex].nodes[1] == node1) &&
           (faceIndices[faceIndex].nodes[2] == node2))
        {
          cellNeighbours[i].faces[0] = faceIndex;
        }

        if((faceIndices[faceIndex].nodes[0] == node0) &&
           (faceIndices[faceIndex].nodes[1] == node1) &&
           (faceIndices[faceIndex].nodes[2] == node3))
        {
          cellNeighbours[i].faces[1] = faceIndex;
        }

        if((faceIndices[faceIndex].nodes[0] == node0) &&
           (faceIndices[faceIndex].nodes[1] == node2) &&
           (faceIndices[faceIndex].nodes[2] == node3))
        {
          cellNeighbours[i].faces[3] = faceIndex;
        }
      }
      for(IndexType j = 0; j < nodeNeighbours[node1].incidentFacesCount; j++)
      {
        IndexType faceIndex = nodeIncidentFacesPool[nodeNeighbours[node1].incidentFacesStart + j];

        if((faceIndices[faceIndex].nodes[0] == node1) &&
           (faceIndices[faceIndex].nodes[1] == node2) &&
           (faceIndices[faceIndex].nodes[2] == node3))
        {
          cellNeighbours[i].faces[2] = faceIndex;
        }
      }
    }
  }
  void BuildCellIncidentFaces()
  {
    for(IndexType i = 0; i < facesCount; i++)
    {
      IndexType faceNode0 = faceIndices[i].nodes[0];
      IndexType faceNode1 = faceIndices[i].nodes[1];
      IndexType faceNode2 = faceIndices[i].nodes[2];
      for(IndexType j = 0; j < nodeNeighbours[faceNode0].incidentCellsCount; j++)
      {
        IndexType cellIndex = nodeIncidentCellsPool[nodeNeighbours[faceNode0].incidentCellsStart + j];

        IndexType cellNode0 = cellIndices[cellIndex].nodes[0];
        IndexType cellNode1 = cellIndices[cellIndex].nodes[1];
        IndexType cellNode2 = cellIndices[cellIndex].nodes[2];
        IndexType cellNode3 = cellIndices[cellIndex].nodes[3];

        IndexType faceInCellIndex = -1;

        if((cellNode0 == faceNode0) && (cellNode1 == faceNode1) && (cellNode2 == faceNode2))
          faceInCellIndex = 0;

        if((cellNode0 == faceNode0) && (cellNode1 == faceNode1) && (cellNode3 == faceNode2))
          faceInCellIndex = 1;

        if((cellNode0 == faceNode0) && (cellNode2 == faceNode1) && (cellNode3 == faceNode2))
          faceInCellIndex = 2;

        if((cellNode1 == faceNode0) && (cellNode2 == faceNode1) && (cellNode3 == faceNode2))
          faceInCellIndex = 3;

        if(faceInCellIndex != IndexType(-1))
          cellNeighbours[cellIndex].faces[faceInCellIndex] = i;
      }
    }
  }

  void BuildEdgeIncidentCells()
  {
    totalEdgeIncidentCells = 0;
    for(IndexType i = 0; i < edgesCount; i++)
    {
      IndexType nodeIndex = edgeIndices[i].nodes[0];

      for(IndexType j = 0; j < nodeNeighbours[nodeIndex].incidentCellsCount; j++)
      {
        IndexType cellIndex = nodeIncidentCellsPool[nodeNeighbours[nodeIndex].incidentCellsStart + j]; 
        for(IndexType k = 0; k < 6; k++)
        {
          if(cellNeighbours[cellIndex].edges[k] == i)
          {
            edgeNeighbours[i].incidentCellsCount++;
            totalEdgeIncidentCells++;
          }
        }
      }
    }

    edgeIncidentCellsPool.resize(totalEdgeIncidentCells);

    IndexType poolState = 0;
    for(IndexType i = 0; i < edgesCount; i++)
    {
      edgeNeighbours[i].incidentCellsStart = poolState;
      poolState += edgeNeighbours[i].incidentCellsCount;
      edgeNeighbours[i].incidentCellsCount = 0;
    }

    if(poolState != totalEdgeIncidentCells)
    {
      //fail
    }

    for(IndexType i = 0; i < edgesCount; i++)
    {
      IndexType nodeIndex = edgeIndices[i].nodes[0];

      for(IndexType j = 0; j < nodeNeighbours[nodeIndex].incidentCellsCount; j++)
      {
        IndexType cellIndex = nodeIncidentCellsPool[nodeNeighbours[nodeIndex].incidentCellsStart + j]; 
        for(IndexType k = 0; k < 6; k++)
        {
          if(cellNeighbours[cellIndex].edges[k] == i)
          {
            edgeIncidentCellsPool[edgeNeighbours[i].incidentCellsStart + edgeNeighbours[i].incidentCellsCount++] = cellIndex;
          }
        }
      }
    }
  }

  void BuildFaceIncidentCells()
  {
    totalFaceIncidentCells = 0;
    for(IndexType i = 0; i < facesCount; i++)
    {
      IndexType nodeIndex = faceIndices[i].nodes[0];

      for(IndexType j = 0; j < nodeNeighbours[nodeIndex].incidentCellsCount; j++)
      {
        IndexType cellIndex = nodeIncidentCellsPool[nodeNeighbours[nodeIndex].incidentCellsStart + j]; 
        for(IndexType k = 0; k < 4; k++)
        {
          if(cellNeighbours[cellIndex].faces[k] == i)
          {
            faceNeighbours[i].incidentCellsCount++;
            totalFaceIncidentCells++;
          }
        }
      }
    }

    faceIncidentCellsPool.resize(totalFaceIncidentCells);

    IndexType poolState = 0;
    for(IndexType i = 0; i < facesCount; i++)
    {
      faceNeighbours[i].incidentCellsStart = poolState;
      poolState += faceNeighbours[i].incidentCellsCount;
      faceNeighbours[i].incidentCellsCount = 0;
    }

    if(poolState != totalFaceIncidentCells)
    {
      // fail
    }

    for(IndexType i = 0; i < facesCount; i++)
    {
      IndexType nodeIndex = faceIndices[i].nodes[0];

      for(IndexType j = 0; j < nodeNeighbours[nodeIndex].incidentCellsCount; j++)
      {
        IndexType cellIndex = nodeIncidentCellsPool[nodeNeighbours[nodeIndex].incidentCellsStart + j]; 
        for(IndexType k = 0; k < 4; k++)
        {
          if(cellNeighbours[cellIndex].faces[k] == i)
          {
            faceIncidentCellsPool[faceNeighbours[i].incidentCellsStart + faceNeighbours[i].incidentCellsCount++] = cellIndex;
          }
        }
      }
    }
  }

public:
  TopologyReconstructor(IndexType *indexData, IndexType indicesCount)
  {
    {
      BuildCellIndices(indexData, indicesCount);
      BuildFaceIndices();
      BuildEdgeIndices();
    }
    {
      InitNeighbouringStructures();
    }
    {
      BuildNodeIncidentCells();
      /*std::cout << "building node incident faces\n";
      BuildNodeIncidentFaces();
      std::cout << "building node incident edges\n";
      BuildNodeIncidentEdges();*/
    }
    {
      BuildCellIncidentEdges();
      BuildCellIncidentFaces();
    }
    {
      BuildEdgeIncidentCells();
      BuildFaceIncidentCells();
    }
    //calculate
  }

  //Nodes

  IndexType GetNodesCount()
  {
    return nodesCount;
  }

  IndexType GetTotalNodeIncidentCells()
  {
    return IndexType(nodeIncidentCellsPool.size());
  }

  IndexType GetNodeIncidentCellsCount(IndexType nodeIndex)
  {
    return nodeNeighbours[nodeIndex].incidentCellsCount;
  }

  void GetNodeIncidentCells(IndexType nodeIndex, IndexType *o_incidentCells)
  {
    for(IndexType i = 0; i < nodeNeighbours[nodeIndex].incidentCellsCount; i++)
    {
      o_incidentCells[i] = nodeIncidentCellsPool[nodeNeighbours[nodeIndex].incidentCellsStart + i];
    }
  }

  //Edges

  IndexType GetEdgesCount()
  {
    return edgesCount;
  }

  void GetEdgeIncidentNodes(IndexType edgeIndex, IndexType *o_incidentNodes)
  {
    o_incidentNodes[0] = edgeIndices[edgeIndex].nodes[0];
    o_incidentNodes[1] = edgeIndices[edgeIndex].nodes[1];
  }

  IndexType GetTotalEdgeIncidentCells()
  {
    return IndexType(edgeIncidentCellsPool.size());
  }

  IndexType GetEdgeIncidentCellsCount(IndexType edgeIndex)
  {
    return edgeNeighbours[edgeIndex].incidentCellsCount;
  }

  void GetEdgeIncidentCells(IndexType edgeIndex, IndexType *o_incidentCells)
  {
    for(IndexType i = 0; i < edgeNeighbours[edgeIndex].incidentCellsCount; i++)
    {
      o_incidentCells[i] = edgeIncidentCellsPool[edgeNeighbours[edgeIndex].incidentCellsStart + i];
    }
  }

  //Faces

  IndexType GetFacesCount()
  {
    return facesCount;
  }

  void GetFaceIncidentNodes(IndexType faceIndex, IndexType *o_incidentNodes)
  {
    o_incidentNodes[0] = faceIndices[faceIndex].nodes[0];
    o_incidentNodes[1] = faceIndices[faceIndex].nodes[1];
    o_incidentNodes[2] = faceIndices[faceIndex].nodes[2];
  }

  IndexType GetTotalFaceIncidentCells()
  {
    return IndexType(faceIncidentCellsPool.size());
  }

  IndexType GetFaceIncidentCellsCount(IndexType faceIndex)
  {
    return faceNeighbours[faceIndex].incidentCellsCount;
  }

  void GetFaceIncidentCells(IndexType faceIndex, IndexType *o_incidentCells)
  {
    for(IndexType i = 0; i < faceNeighbours[faceIndex].incidentCellsCount; i++)
    {
      o_incidentCells[i] = faceIncidentCellsPool[faceNeighbours[faceIndex].incidentCellsStart + i];
    }
  }

  //Cells

  IndexType GetCellsCount()
  {
    return cellsCount;
  }

  void GetCellIncidentNodes(IndexType cellIndex, IndexType *o_nodes)
  {
    o_nodes[0] = cellIndices[cellIndex].nodes[0];
    o_nodes[1] = cellIndices[cellIndex].nodes[1];
    o_nodes[2] = cellIndices[cellIndex].nodes[2];
    o_nodes[3] = cellIndices[cellIndex].nodes[3];
  }

  void GetCellIncidentEdges(IndexType cellIndex, IndexType *o_edges)
  {
    o_edges[0] = cellNeighbours[cellIndex].edges[0];
    o_edges[1] = cellNeighbours[cellIndex].edges[1];
    o_edges[2] = cellNeighbours[cellIndex].edges[2];
    o_edges[3] = cellNeighbours[cellIndex].edges[3];
    o_edges[4] = cellNeighbours[cellIndex].edges[4];
    o_edges[5] = cellNeighbours[cellIndex].edges[5];
  }

  void GetCellIncidentFaces(IndexType cellIndex, IndexType *o_faces)
  {
    o_faces[0] = cellNeighbours[cellIndex].faces[0];
    o_faces[1] = cellNeighbours[cellIndex].faces[1];
    o_faces[2] = cellNeighbours[cellIndex].faces[2];
    o_faces[3] = cellNeighbours[cellIndex].faces[3];
  }
private:
  std::vector<NodeNeighbours> nodeNeighbours;
  std::vector<EdgeNeighbours> edgeNeighbours;
  std::vector<FaceNeighbours> faceNeighbours;
  std::vector<CellNeighbours> cellNeighbours;

  std::vector<CellIndices> cellIndices;

  DuplicateRemover<EdgeIndices, IndexType> edgeIndices;
  DuplicateRemover<FaceIndices, IndexType> faceIndices;

  std::vector<IndexType> nodeIncidentCellsPool;
  IndexType totalNodeIncidentCells;
  std::vector<IndexType> nodeIncidentFacesPool;
  IndexType totalNodeIncidentFaces;
  std::vector<IndexType> nodeIncidentEdgesPool;
  IndexType totalNodeIncidentEdges;

  std::vector<IndexType> edgeIncidentCellsPool;
  IndexType totalEdgeIncidentCells;

  std::vector<IndexType> faceIncidentCellsPool;
  IndexType totalFaceIncidentCells;

  IndexType cellsCount;
  IndexType facesCount;
  IndexType edgesCount;
  IndexType nodesCount;
};

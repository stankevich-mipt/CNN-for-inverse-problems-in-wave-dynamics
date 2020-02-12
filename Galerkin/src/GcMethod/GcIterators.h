template <typename ElasticMesh>
class NodePointGcIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;
  IndexType pointIndex;
  ElasticMesh *mesh;

  NodePointGcIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  NodePointGcIterator<ElasticMesh> operator + (IndexType offset)
  {
    return NodePointGcIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (NodePointGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
    return mesh->GetNodeElasticData(pointIndex)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return mesh->GetNodeElasticData(pointIndex)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    return mesh->GetNodePointPos(pointIndex);
//    return reconstructor->GetMesh()->GetNode(pointIndex)->pos;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    posCell = mesh->GetCellInNodeDirection(pointIndex,  dir);
    negCell = mesh->GetCellInNodeDirection(pointIndex, -dir);
  }
};


template <typename ElasticMesh>
class EdgePointGcIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType pointsCount = ElasticMesh::edgePointsCount;

  EdgePointGcIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  EdgePointGcIterator<ElasticMesh> operator + (IndexType offset)
  {
    return EdgePointGcIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (EdgePointGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
    return mesh->GetEdgeElasticData(pointIndex)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return mesh->GetEdgeElasticData(pointIndex)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    return mesh->GetEdgePointPos(pointIndex / pointsCount, pointIndex % pointsCount);
//    return reconstructor->GetMesh()->GetEdge(pointIndex / edgePointsCount)->data.controlPoints[pointIndex % edgePointsCount].pos;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    posCell = mesh->GetCellInEdgeDirection(pointIndex / pointsCount,  dir);
    negCell = mesh->GetCellInEdgeDirection(pointIndex / pointsCount, -dir);
  }
};

template <typename ElasticMesh>
class FacePointGcIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType pointsCount = ElasticMesh::facePointsCount;

  FacePointGcIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  FacePointGcIterator<ElasticMesh> operator + (IndexType offset)
  {
    return FacePointGcIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (FacePointGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
    return mesh->GetFaceElasticData(pointIndex)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return mesh->GetFaceElasticData(pointIndex)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    return mesh->GetFacePointPos(pointIndex / pointsCount, pointIndex % pointsCount);
//    return reconstructor->GetMesh()->GetEdge(pointIndex / edgePointsCount)->data.controlPoints[pointIndex % edgePointsCount].pos;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    if(pointsCount > 0)
    {
      posCell = mesh->GetCellInFaceDirection(pointIndex / pointsCount,  dir);
      negCell = mesh->GetCellInFaceDirection(pointIndex / pointsCount, -dir);
    }else
    {
      posCell = -1;
      negCell = -1;
    }
  }
};

template <typename ElasticMesh>
class CellPointGcIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType pointsCount = ElasticMesh::cellPointsCount;

  CellPointGcIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  CellPointGcIterator<ElasticMesh> operator + (IndexType offset)
  {
    return CellPointGcIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (CellPointGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
    return mesh->GetCellElasticData(pointIndex)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return mesh->GetCellElasticData(pointIndex)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    return mesh->GetCellPointPos(pointIndex / pointsCount, pointIndex % pointsCount);
//    return reconstructor->GetMesh()->GetEdge(pointIndex / edgePointsCount)->data.controlPoints[pointIndex % edgePointsCount].pos;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    posCell = pointIndex / pointsCount; //pointsCount per each cell
    negCell = pointIndex / pointsCount;
  }
};

template <typename ElasticMesh>
class UniversalGcIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  IndexType submeshIndex;

  static const IndexType cellPointsCount = ElasticMesh::cellPointsCount;
  static const IndexType facePointsCount = ElasticMesh::facePointsCount;
  static const IndexType edgePointsCount = ElasticMesh::edgePointsCount;

  UniversalGcIterator(ElasticMesh *i_mesh, IndexType i_pointNumber, IndexType i_submeshIndex)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
    submeshIndex = i_submeshIndex;
  }
public:
  UniversalGcIterator<ElasticMesh> operator + (IndexType offset)
  {
    return UniversalGcIterator<ElasticMesh>(mesh, pointIndex + offset, submeshIndex);
  }

  IndexType operator - (UniversalGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  IndexType GetSubmeshIndex()
  {
    return submeshIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {

    return GetElasticData()->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return GetElasticData()->nextLayer;
  }

  typename ElasticMesh::Elastic &tempValue()
  {
    return GetElasticData()->tempLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        return mesh->GetSubmeshNodePointPos(pointIndex - GetNodeOffset(), submeshIndex);
      }break;

      case Edge:
      {
        return mesh->GetSubmeshEdgePointPos(pointIndex - GetEdgeOffset(), submeshIndex);
      }break;

      case Face:
      {
        return mesh->GetSubmeshFacePointPos(pointIndex - GetFaceOffset(), submeshIndex);
      }break;

      case Cell:
      {
        return mesh->GetSubmeshCellPointPos(pointIndex - GetCellOffset(), submeshIndex);
      }break;
    }
    return typename ElasticMesh::Vector3(-1.0, -1.0, -1.0);
  }

  typename ElasticMesh::ElasticPoint *GetElasticData()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        return mesh->GetSubmeshNodeElasticData(pointIndex - GetNodeOffset(), submeshIndex);
      }break;

      case Edge:
      {
        return mesh->GetSubmeshEdgeElasticData(pointIndex - GetEdgeOffset(), submeshIndex);
      }break;

      case Face:
      {
        return mesh->GetSubmeshFaceElasticData(pointIndex - GetFaceOffset(), submeshIndex);
      }break;

      case Cell:
      {
        return mesh->GetSubmeshCellElasticData(pointIndex - GetCellOffset(), submeshIndex);
      }break;
    }
    return 0;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeIndex = mesh->GetSubmeshNodeGlobalIndex(pointIndex - GetNodeOffset(), submeshIndex);
        posCell = mesh->GetCellInNodeDirection(nodeIndex, dir);
        negCell = mesh->GetCellInNodeDirection(nodeIndex, -dir);
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgeIndex = mesh->GetSubmeshEdgeGlobalIndex((pointIndex - GetEdgeOffset()) / edgePointsCount, submeshIndex);
          posCell = mesh->GetCellInEdgeDirection(edgeIndex,  dir);
          negCell = mesh->GetCellInEdgeDirection(edgeIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType faceIndex = mesh->GetSubmeshFaceGlobalIndex((pointIndex - GetFaceOffset()) / facePointsCount, submeshIndex);
          posCell = mesh->GetCellInFaceDirection(faceIndex,  dir);
          negCell = mesh->GetCellInFaceDirection(faceIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;

      case Cell:
      {
        if(cellPointsCount > 0)
        {
          IndexType cellIndex = mesh->GetSubmeshCellGlobalIndex((pointIndex - GetCellOffset()) / cellPointsCount, submeshIndex);
          posCell = cellIndex;
          negCell = cellIndex;
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;
    }
  }
private:
  enum FeatureType
  {
    Node,
    Edge,
    Face,
    Cell
  };
  IndexType GetNodeOffset()
  {
    return 0;
  }
  IndexType GetEdgeOffset()
  {
    return GetNodeOffset() + mesh->GetSubmeshNodesCount(submeshIndex) * 1;
  }
  IndexType GetFaceOffset()
  {
    return GetEdgeOffset() + mesh->GetSubmeshEdgesCount(submeshIndex) * edgePointsCount;
  }
  IndexType GetCellOffset()
  {
    return GetFaceOffset() + mesh->GetSubmeshFacesCount(submeshIndex) * facePointsCount;
  }

  FeatureType GetFeatureTypeByPointIndex(IndexType pointIndex)
  {
    if(pointIndex < GetEdgeOffset()) return Node;
    if(pointIndex < GetFaceOffset()) return Edge;
    if(pointIndex < GetCellOffset()) return Face;
    return Cell;
  }
};




/*template <typename ElasticMesh>
class NodePointBorderIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  NodePointBorderIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  NodePointBorderIterator<ElasticMesh> operator + (IndexType offset)
  {
    return NodePointBorderIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (NodePointBorderIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
    return mesh->GetBorderNodeElasticData(pointIndex)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return mesh->GetBorderNodeElasticData(pointIndex)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    return mesh->GetBorderNodePointPos(pointIndex);
//    return reconstructor->GetMesh()->GetNode(pointIndex)->pos;
  }

  typename ElasticMesh::Vector3 currNormal()
  {
    return mesh->GetBorderNodeNormal(pointIndex);
//    return reconstructor->GetMesh()->GetNode(pointIndex)->pos;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    posCell = mesh->GetCellInNodeDirection(mesh->GetBorderNodeIndex(pointIndex),  dir);
    negCell = mesh->GetCellInNodeDirection(mesh->GetBorderNodeIndex(pointIndex), -dir);
  }
};

template <typename ElasticMesh>
class EdgePointBorderIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType pointsCount = ElasticMesh::edgePointsCount;

  EdgePointBorderIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  EdgePointBorderIterator<ElasticMesh> operator + (IndexType offset)
  {
    return EdgePointBorderIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (EdgePointBorderIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
//    return mesh->GetEdgeElasticData(mesh->GetBorderEdgeIndex(pointIndex / pointsCount), pointIndex % pointsCount)->currLayer;
    return mesh->GetBorderEdgeElasticData(pointIndex / pointsCount, pointIndex % pointsCount)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return mesh->GetBorderEdgeElasticData(pointIndex / pointsCount, pointIndex % pointsCount)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    return mesh->GetEdgePointPos(mesh->GetBorderEdgeIndex(pointIndex / pointsCount), pointIndex % pointsCount);
//    return reconstructor->GetMesh()->GetEdge(pointIndex / edgePointsCount)->data.controlPoints[pointIndex % edgePointsCount].pos;
  }

  typename ElasticMesh::Vector3 currNormal()
  {
    return mesh->GetBorderEdgePointNormal(pointIndex / pointsCount, pointIndex % pointsCount);
//    return reconstructor->GetMesh()->GetEdge(pointIndex / edgePointsCount)->data.controlPoints[pointIndex % edgePointsCount].pos;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    posCell = mesh->GetCellInEdgeDirection(mesh->GetBorderEdgeIndex(pointIndex / pointsCount),  dir);
    negCell = mesh->GetCellInEdgeDirection(mesh->GetBorderEdgeIndex(pointIndex / pointsCount), -dir);
  }
};


template <typename ElasticMesh>
class FacePointBorderIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType pointsCount = ElasticMesh::facePointsCount;

  FacePointBorderIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  FacePointBorderIterator<ElasticMesh> operator + (IndexType offset)
  {
    return FacePointBorderIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (FacePointBorderIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
    return mesh->GetBorderFaceElasticData(pointIndex / pointsCount, pointIndex % pointsCount)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return mesh->GetBorderFaceElasticData(pointIndex / pointsCount, pointIndex % pointsCount)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    return mesh->GetBorderFacePointPos(pointIndex / pointsCount, pointIndex % pointsCount);
  }

  typename ElasticMesh::Vector3 currNormal()
  {
    return mesh->GetBorderFacePointNormal(pointIndex / pointsCount, pointIndex % pointsCount);
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    posCell = mesh->GetCellInFaceDirection(mesh->GetBorderFaceIndex(pointIndex / pointsCount),  dir);
    negCell = mesh->GetCellInFaceDirection(mesh->GetBorderFaceIndex(pointIndex / pointsCount), -dir);
  }
};


template <typename ElasticMesh>
class NodePointContactIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  NodePointContactIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  NodePointContactIterator<ElasticMesh> operator + (IndexType offset)
  {
    return NodePointContactIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (NodePointContactIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue0()
  {
    return mesh->GetContactNodeElasticPair(pointIndex).elasticPoints[0]->currLayer;
  }

  typename ElasticMesh::Elastic &currValue1()
  {
    return mesh->GetContactNodeElasticPair(pointIndex).elasticPoints[1]->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue0()
  {
    return mesh->GetContactNodeElasticPair(pointIndex).elasticPoints[0]->nextLayer;
  }

  typename ElasticMesh::Elastic &nextValue1()
  {
    return mesh->GetContactNodeElasticPair(pointIndex).elasticPoints[1]->nextLayer;
  }

  typename ElasticMesh::Vector3 currNormal()
  {
    return mesh->GetContactNodeNormals(pointIndex).normals[0];
  }
};

template <typename ElasticMesh>
class EdgePointContactIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType pointsCount = ElasticMesh::edgePointsCount;

  EdgePointContactIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  EdgePointContactIterator<ElasticMesh> operator + (IndexType offset)
  {
    return EdgePointContactIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (EdgePointContactIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue0()
  {
    return mesh->GetContactEdgeElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[0]->currLayer;
  }

  typename ElasticMesh::Elastic &currValue1()
  {
    return mesh->GetContactEdgeElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[1]->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue0()
  {
    return mesh->GetContactEdgeElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[0]->nextLayer;
  }

  typename ElasticMesh::Elastic &nextValue1()
  {
    return mesh->GetContactEdgeElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[1]->nextLayer;
  }


  typename ElasticMesh::Vector3 currNormal()
  {
    return mesh->GetContactEdgePointNormal(pointIndex / pointsCount, pointIndex % pointsCount);
//    return reconstructor->GetMesh()->GetEdge(pointIndex / edgePointsCount)->data.controlPoints[pointIndex % edgePointsCount].pos;
  }

};

template <typename ElasticMesh>
class FacePointContactIterator
{
public:
  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType pointsCount = ElasticMesh::facePointsCount;

  FacePointContactIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  FacePointContactIterator<ElasticMesh> operator + (IndexType offset)
  {
    return FacePointContactIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (FacePointContactIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue0()
  {
    return mesh->GetContactFaceElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[0]->currLayer;
  }

  typename ElasticMesh::Elastic &currValue1()
  {
    return mesh->GetContactFaceElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[1]->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue0()
  {
    return mesh->GetContactFaceElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[0]->nextLayer;
  }

  typename ElasticMesh::Elastic &nextValue1()
  {
    return mesh->GetContactFaceElasticPair(pointIndex / pointsCount, pointIndex % pointsCount).elasticPoints[1]->nextLayer;
  }


  typename ElasticMesh::Vector3 currNormal()
  {
    return mesh->GetContactFacePointNormal(pointIndex / pointsCount, pointIndex % pointsCount);
  }

};

template <typename ElasticMesh>
class UniversalGcBorderIterator
{
public:
  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType cellPointsCount = ElasticMesh::cellPointsCount;
  static const IndexType facePointsCount = ElasticMesh::facePointsCount;
  static const IndexType edgePointsCount = ElasticMesh::edgePointsCount;

  UniversalGcBorderIterator(){}
  UniversalGcBorderIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  UniversalGcBorderIterator<ElasticMesh> operator + (IndexType offset)
  {
    return UniversalGcBorderIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (UniversalGcBorderIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {

    return GetElasticData()->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return GetElasticData()->nextLayer;
  }

  typename ElasticMesh::ElasticPoint *GetElasticData()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetBorderNodeElasticData(nodePointIndex);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetBorderEdgeElasticData(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetBorderFaceElasticData(facePointIndex / facePointsCount, facePointIndex % facePointsCount);
      }break;
    }
    assert(0);
    return 0;
  }

  typename ElasticMesh::Vector3 currNormal()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetBorderNodeNormal(nodePointIndex);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetBorderEdgePointNormal(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetBorderFacePointNormal(facePointIndex / facePointsCount, facePointIndex % facePointsCount);
      }break;
    }
    assert(0);
    return ElasticMesh::Vector3(0, 0, 0);
  }

  IndexType GetSubmeshIndex()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeIndex = mesh->GetBorderNodeIndex(pointIndex - GetNodeOffset());
        return mesh->GetNodeSubmeshIndex(nodeIndex);
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgeIndex = mesh->GetBorderEdgeIndex((pointIndex - GetEdgeOffset()) / edgePointsCount);
          return mesh->GetEdgeSubmeshIndex(edgeIndex);
        }else
        {
          assert(0);
        }
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType faceIndex = mesh->GetBorderFaceIndex((pointIndex - GetFaceOffset()) / facePointsCount);
          return mesh->GetFaceSubmeshIndex(faceIndex);
        }else
        {
          assert(0);
        }
      }break;
    }
    assert(0);
    return -1;
  }

private:
  enum FeatureType
  {
    Node,
    Edge,
    Face,
    Cell
  };
  IndexType GetNodeOffset()
  {
    return 0;
  }
  IndexType GetEdgeOffset()
  {
    return GetNodeOffset() + mesh->GetBorderNodesCount() * 1;
  }
  IndexType GetFaceOffset()
  {
    return GetEdgeOffset() + mesh->GetBorderEdgesCount() * edgePointsCount;
  }


  FeatureType GetFeatureTypeByPointIndex(IndexType pointIndex)
  {
    if(pointIndex < GetEdgeOffset()) return Node;
    if(pointIndex < GetFaceOffset()) return Edge;
    return Face;
  }
};

template <typename ElasticMesh>
class UniversalGcContactIterator
{
public:
  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType cellPointsCount = ElasticMesh::cellPointsCount;
  static const IndexType facePointsCount = ElasticMesh::facePointsCount;
  static const IndexType edgePointsCount = ElasticMesh::edgePointsCount;

  UniversalGcContactIterator(){}
  UniversalGcContactIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  UniversalGcContactIterator<ElasticMesh> operator + (IndexType offset)
  {
    return UniversalGcContactIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (UniversalGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue(IndexType contactSide)
  {

    return GetElasticData(contactSide)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue(IndexType contactSide)
  {
    return GetElasticData(contactSide)->nextLayer;
  }

  typename ElasticMesh::Vector3 currPos(IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodePositionsPair(nodePointIndex).vectors[contactSide];
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetContactEdgePositionsPair(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount).vectors[contactSide];
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetContactFacePositionsPair(facePointIndex / facePointsCount, facePointIndex % facePointsCount).vectors[contactSide];
      }break;
    }
    assert(0);
    return ElasticMesh::Vector3(0, 0, 0);
  }

  typename ElasticMesh::ElasticPoint *GetElasticData(IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodeElasticPair(nodePointIndex).elasticPoints[contactSide];
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetContactEdgeElasticPair(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount).elasticPoints[contactSide];
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetContactFaceElasticPair(facePointIndex / facePointsCount, facePointIndex % facePointsCount).elasticPoints[contactSide];
      }break;
    }
    assert(0);
    return 0;
  }

  typename ElasticMesh::Vector3 currNormal()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodeNormals(nodePointIndex).normals[0];
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetContactEdgePointNormal(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetContactFacePointNormal(facePointIndex / facePointsCount, facePointIndex % facePointsCount);
      }break;
    }
    assert(0);
    return ElasticMesh::Vector3(0, 0, 0);
  }

  IndexType GetSubmeshIndex(IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeIndex = mesh->GetContactNodePair(pointIndex - GetNodeOffset()).nodes[contactSide];
        return mesh->GetNodeSubmeshIndex(nodeIndex);
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgeIndex = mesh->GetContactEdgePair((pointIndex - GetEdgeOffset()) / edgePointsCount).edges[contactSide];
          return mesh->GetEdgeSubmeshIndex(edgeIndex);
        }else
        {
          assert(0);
        }
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType faceIndex = mesh->GetContactFacePair((pointIndex - GetFaceOffset()) / facePointsCount).faces[contactSide];
          return mesh->GetFaceSubmeshIndex(faceIndex);
        }else
        {
          assert(0);
        }
      }break;
    }
    assert(0);
    return -1;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell, IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeIndex = mesh->GetContactNodePair(pointIndex - GetNodeOffset()).nodes[contactSide];
        posCell = mesh->GetCellInNodeDirection(nodeIndex, dir);
        negCell = mesh->GetCellInNodeDirection(nodeIndex, -dir);
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgeIndex = mesh->GetContactEdgePair((pointIndex - GetEdgeOffset()) / edgePointsCount).edges[contactSide];
          posCell = mesh->GetCellInEdgeDirection(edgeIndex,  dir);
          negCell = mesh->GetCellInEdgeDirection(edgeIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType faceIndex = mesh->GetContactFacePair((pointIndex - GetFaceOffset()) / facePointsCount).faces[contactSide];
          posCell = mesh->GetCellInFaceDirection(faceIndex,  dir);
          negCell = mesh->GetCellInFaceDirection(faceIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;
    }
    assert(0);
  }
  struct SideContactPointGcIterator
  {
    SideContactPointGcIterator(UniversalGcContactIterator<ElasticMesh> contactIterator, IndexType contactSide)
    {
      this->contactIterator = contactIterator;
      this->contactSide = contactSide;
      tmpPoint = nextValue();
    }
    SideContactPointGcIterator(){}
    SideContactPointGcIterator operator + (IndexType offset)
    {
      return SideContactPointGcIterator(this->contactIterator + offset, contactSide);
    }

    IndexType operator - (SideContactPointGcIterator &it)
    {
      return 1;//this->contactIterator - it.contactIterator;
    }
    void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
    {
      this->contactIterator.findCellAlongDir(dir, negCell, posCell, contactSide);
    }

    typename ElasticMesh::Elastic &currValue()
    {

      return this->contactIterator.currValue(contactSide);
    }

    typename ElasticMesh::Elastic &nextValue()
    {
      return tmpPoint;
    }

    typename ElasticMesh::Vector3 currPos()
    {
      return this->contactIterator.currPos(contactSide);
    }
  private:
    typename ElasticMesh::Elastic tmpPoint;
    UniversalGcContactIterator<ElasticMesh> contactIterator;
    IndexType contactSide;
  };

  SideContactPointGcIterator GetSidePointIterator(IndexType contactSide)
  {
    return SideContactPointGcIterator(*this, contactSide);
  }
private:
  enum FeatureType
  {
    Node,
    Edge,
    Face,
    Cell
  };
  IndexType GetNodeOffset()
  {
    return 0;
  }
  IndexType GetEdgeOffset()
  {
    return GetNodeOffset() + mesh->GetContactNodesCount() * 1;
  }
  IndexType GetFaceOffset()
  {
    return GetEdgeOffset() + mesh->GetContactEdgesCount() * edgePointsCount;
  }


  FeatureType GetFeatureTypeByPointIndex(IndexType pointIndex)
  {
    if(pointIndex < GetEdgeOffset()) return Node;
    if(pointIndex < GetFaceOffset()) return Edge;
    return Face;
  }
};*/


template <typename ElasticMesh>
class ContactGroupGcIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;

  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType cellPointsCount = ElasticMesh::cellPointsCount;
  static const IndexType facePointsCount = ElasticMesh::facePointsCount;
  static const IndexType edgePointsCount = ElasticMesh::edgePointsCount;

  ContactGroupGcIterator(){}
  ContactGroupGcIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  ContactGroupGcIterator<ElasticMesh> operator + (IndexType offset)
  {
    return ContactGroupGcIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (ContactGroupGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue(IndexType contactSide)
  {

    return GetElasticData(contactSide)->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue(IndexType contactSide)
  {
    return GetElasticData(contactSide)->nextLayer;
  }

  typename ElasticMesh::Elastic &tempValue(IndexType contactSide)
  {
    return GetElasticData(contactSide)->tempLayer;
  }

  typename ElasticMesh::Vector3 currPos(IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeGroupIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodeGroupPointPos(nodeGroupIndex, contactSide);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetContactEdgeGroupPointPos(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount, contactSide);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetContactFaceGroupPointPos(facePointIndex / facePointsCount, facePointIndex % facePointsCount, contactSide);
      }break;
    }
    assert(0);
    return typename ElasticMesh::Vector3(0, 0, 0);
  }

  typename ElasticMesh::ElasticPoint *GetElasticData(IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodeGroupElasticData(nodePointIndex, contactSide);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetContactEdgeGroupElasticData(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount, contactSide);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetContactFaceGroupElasticData(facePointIndex / facePointsCount, facePointIndex % facePointsCount, contactSide);
      }break;
    }
    assert(0);
    return 0;
  }

  typename ElasticMesh::Vector3 currNormal(IndexType contactSide0, IndexType contactSide1)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodePairInfo(nodePointIndex, contactSide0, contactSide1).normal;
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetContactEdgePointPairInfo(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount, contactSide0, contactSide1).normal;
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetContactFacePointPairInfo(facePointIndex / facePointsCount, facePointIndex % facePointsCount, contactSide0, contactSide1).normal;
      }break;
    }
    assert(0);
    return typename ElasticMesh::Vector3(0, 0, 0);
  }

  IndexType GetContactTypeIndex(IndexType contactSide0, IndexType contactSide1)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodePairInfo(nodePointIndex, contactSide0, contactSide1).typeIndex;
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetContactEdgePointPairInfo(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount, contactSide0, contactSide1).typeIndex;
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetContactFacePointPairInfo(facePointIndex / facePointsCount, facePointIndex % facePointsCount, contactSide0, contactSide1).typeIndex;
      }break;
    }
    assert(0);
    return IndexType(-1);
  }

  IndexType GetSubmeshIndex(IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetNodeSubmeshIndex(mesh->GetContactNodeInGroup(nodePointIndex, contactSide).nodeIndex);
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgePointIndex = pointIndex - GetEdgeOffset();
          return mesh->GetEdgeSubmeshIndex(mesh->GetContactEdgeInGroup(edgePointIndex / edgePointsCount, contactSide).edgeIndex);
        }else
        {
          assert(0);
        }
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType facePointIndex = pointIndex - GetFaceOffset();
          return mesh->GetFaceSubmeshIndex(mesh->GetContactFaceInGroup(facePointIndex / facePointsCount, contactSide).faceIndex);
        }else
        {
          assert(0);
        }
      }break;
    }
    assert(0);
    return -1;
  }

  IndexType GetContactSidesCount()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetContactNodeGroupSize(nodePointIndex);
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgePointIndex = pointIndex - GetEdgeOffset();
          return mesh->GetContactEdgeGroupSize(edgePointIndex / edgePointsCount);
        }else
        {
          assert(0);
        }
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType facePointIndex = pointIndex - GetFaceOffset();
          return mesh->GetContactFaceGroupSize(facePointIndex / facePointsCount);
        }else
        {
          assert(0);
        }
      }break;
    }
    assert(0);
    return -1;
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell, IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeIndex = mesh->GetContactNodeInGroup(pointIndex - GetNodeOffset(), contactSide).nodeIndex;
        posCell = mesh->GetCellInNodeDirection(nodeIndex, dir);
        negCell = mesh->GetCellInNodeDirection(nodeIndex, -dir);
        return;
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgeIndex = mesh->GetContactEdgeInGroup((pointIndex - GetEdgeOffset()) / edgePointsCount, contactSide).edgeIndex;
          posCell = mesh->GetCellInEdgeDirection(edgeIndex,  dir);
          negCell = mesh->GetCellInEdgeDirection(edgeIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType faceIndex = mesh->GetContactFaceInGroup((pointIndex - GetFaceOffset()) / facePointsCount, contactSide).faceIndex;
          posCell = mesh->GetCellInFaceDirection(faceIndex,  dir);
          negCell = mesh->GetCellInFaceDirection(faceIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;
    }
    assert(0);
  }

  struct SideContactPointGcIterator
  {
    SideContactPointGcIterator(ContactGroupGcIterator<ElasticMesh> contactIterator, IndexType contactSide)
    {
      this->contactIterator = contactIterator;
      this->contactSide = contactSide;
      tmpPoint = nextValue();
    }
    SideContactPointGcIterator(){}
    SideContactPointGcIterator operator + (IndexType offset)
    {
      return SideContactPointGcIterator(this->contactIterator + offset, contactSide);
    }

    IndexType operator - (SideContactPointGcIterator &it)
    {
      return 1;//this->contactIterator - it.contactIterator;
    }
    void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
    {
      this->contactIterator.findCellAlongDir(dir, negCell, posCell, contactSide);
    }

    typename ElasticMesh::Elastic &currValue()
    {

      return this->contactIterator.currValue(contactSide);
    }

    typename ElasticMesh::Elastic &nextValue()
    {
      return tmpPoint;
    }

    typename ElasticMesh::Vector3 currPos()
    {
      return this->contactIterator.currPos(contactSide);
    }
  private:
    typename ElasticMesh::Elastic tmpPoint;
    ContactGroupGcIterator<ElasticMesh> contactIterator;
    IndexType contactSide;
  };

  SideContactPointGcIterator GetSidePointIterator(IndexType contactSide)
  {
    return SideContactPointGcIterator(*this, contactSide);
  }
private:
  enum FeatureType
  {
    Node,
    Edge,
    Face,
    Cell
  };
  IndexType GetNodeOffset()
  {
    return 0;
  }
  IndexType GetEdgeOffset()
  {
    return GetNodeOffset() + mesh->GetContactNodeGroupsCount() * 1;
  }
  IndexType GetFaceOffset()
  {
    return GetEdgeOffset() + mesh->GetContactEdgeGroupsCount() * edgePointsCount;
  }


  FeatureType GetFeatureTypeByPointIndex(IndexType pointIndex)
  {
    if(pointIndex < GetEdgeOffset()) return Node;
    if(pointIndex < GetFaceOffset()) return Edge;
    return Face;
  }
};


template <typename ElasticMesh>
class GlobalGcIterator
{
public:
  typedef typename ElasticMesh::IndexType IndexType;
  IndexType pointIndex;
  ElasticMesh *mesh;

  static const IndexType cellPointsCount = ElasticMesh::cellPointsCount;
  static const IndexType facePointsCount = ElasticMesh::facePointsCount;
  static const IndexType edgePointsCount = ElasticMesh::edgePointsCount;

  GlobalGcIterator(){}
  GlobalGcIterator(ElasticMesh *i_mesh, IndexType i_pointNumber)
  {
    mesh   = i_mesh;
    pointIndex      = i_pointNumber;
  }
public:
  GlobalGcIterator<ElasticMesh> operator + (IndexType offset)
  {
    return GlobalGcIterator<ElasticMesh>(mesh, pointIndex + offset);
  }

  IndexType operator - (GlobalGcIterator<ElasticMesh> &it)
  {
    return pointIndex - it.pointIndex;
  }

  typename ElasticMesh::Elastic &currValue()
  {
    return GetElasticData()->currLayer;
  }

  typename ElasticMesh::Elastic &nextValue()
  {
    return GetElasticData()->nextLayer;
  }

  typename ElasticMesh::Elastic &tempValue()
  {
    return GetElasticData()->tempLayer;
  }

  typename ElasticMesh::Vector3 currPos()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeIndex = pointIndex - GetNodeOffset();
        return mesh->GetNodePointPos(nodeIndex);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetEdgePointPos(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetFacePointPos(facePointIndex / facePointsCount, facePointIndex % facePointsCount);
      }break;

      case Cell:
      {
        IndexType cellPointIndex = pointIndex - GetCellOffset();
        return mesh->GetCellPointPos(cellPointIndex / cellPointsCount, cellPointIndex % cellPointsCount);
      }break;
    }
    assert(0);
    return typename ElasticMesh::Vector3(0, 0, 0);
  }

  typename ElasticMesh::ElasticPoint *GetElasticData()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetNodeElasticData(nodePointIndex);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetEdgeElasticData(edgePointIndex / edgePointsCount, edgePointIndex % edgePointsCount);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetFaceElasticData(facePointIndex / facePointsCount, facePointIndex % facePointsCount);
      }break;

      case Cell:
      {
        IndexType cellPointIndex = pointIndex - GetCellOffset();
        return mesh->GetCellElasticData(cellPointIndex / cellPointsCount, cellPointIndex % cellPointsCount);
      }break;
     }
    assert(0);
    return 0;
  }

  IndexType GetSubmeshIndex()
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetNodeSubmeshIndex(nodePointIndex);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetEdgeSubmeshIndex(edgePointIndex / edgePointsCount);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetFaceSubmeshIndex(facePointIndex / facePointsCount);
      }break;

      case Cell:
      {
        IndexType cellPointIndex = pointIndex - GetCellOffset();
        return mesh->GetCellSubmeshIndex(cellPointIndex / cellPointsCount);
      }break;
    }
    assert(0);
    return 0;
  }


  typename ElasticMesh::Vector3 currNormal(IndexType contactSide)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodePointIndex = pointIndex - GetNodeOffset();
        return mesh->GetNodeSubmeshIndex(nodePointIndex);
      }break;

      case Edge:
      {
        IndexType edgePointIndex = pointIndex - GetEdgeOffset();
        return mesh->GetEdgeSubmeshIndex(edgePointIndex);
      }break;

      case Face:
      {
        IndexType facePointIndex = pointIndex - GetFaceOffset();
        return mesh->GetFaceSubmeshIndex(facePointIndex);
      }break;

      case Cell:
      {
        IndexType cellPointIndex = pointIndex - GetCellOffset();
        return mesh->GetCellSubmeshIndex(cellPointIndex);
      }break;
    }
    assert(0);
    return typename ElasticMesh::Vector3(0, 0, 0);
  }

  void findCellAlongDir(const typename ElasticMesh::Vector3 &dir, IndexType &negCell, IndexType &posCell)
  {
    switch(GetFeatureTypeByPointIndex(pointIndex))
    {
      case Node:
      {
        IndexType nodeIndex = pointIndex - GetNodeOffset();
        posCell = mesh->GetCellInNodeDirection(nodeIndex, dir);
        negCell = mesh->GetCellInNodeDirection(nodeIndex, -dir);
        return;
      }break;

      case Edge:
      {
        if(edgePointsCount > 0)
        {
          IndexType edgeIndex = (pointIndex - GetEdgeOffset()) / edgePointsCount;
          posCell = mesh->GetCellInEdgeDirection(edgeIndex,  dir);
          negCell = mesh->GetCellInEdgeDirection(edgeIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;

      case Face:
      {
        if(facePointsCount > 0)
        {
          IndexType faceIndex = (pointIndex - GetFaceOffset()) / facePointsCount;
          posCell = mesh->GetCellInFaceDirection(faceIndex,  dir);
          negCell = mesh->GetCellInFaceDirection(faceIndex, -dir);
        }else
        {
          posCell = -1;
          negCell = -1;
        }
        return;
      }break;

      case Cell:
      {
        if(cellPointsCount > 0)
        {
          IndexType cellIndex = (pointIndex - GetCellOffset()) / cellPointsCount;
          posCell = cellIndex;
          negCell = cellIndex;
        }
        return;
      }break;
    }
    assert(0);
  }

private:
  enum FeatureType
  {
    Node,
    Edge,
    Face,
    Cell
  };
  IndexType GetNodeOffset()
  {
    return 0;
  }

  IndexType GetEdgeOffset()
  {
    return GetNodeOffset() + mesh->GetNodesCount() * 1;
  }

  IndexType GetFaceOffset()
  {
    return GetEdgeOffset() + mesh->GetEdgesCount() * edgePointsCount;
  }

  IndexType GetCellOffset()
  {
    return GetFaceOffset() + mesh->GetFacesCount() * facePointsCount;
  }


  FeatureType GetFeatureTypeByPointIndex(IndexType pointIndex)
  {
    if(pointIndex < GetEdgeOffset()) return Node;
    if(pointIndex < GetFaceOffset()) return Edge;
    if(pointIndex < GetCellOffset()) return Face;
    return Cell;
  }
};

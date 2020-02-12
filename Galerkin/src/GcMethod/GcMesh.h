#include "Interpolator.h"
#include "GcIterators.h"
//#include "tetrageom/src/tetrageom.h"
//#include "pointReconstructor.h"


template<typename Space, int splitScheme>
struct ElasticPointT;

template<typename Space>
struct ElasticPointT<Space, 0>
{
  static const typename Space::IndexType paramCount = 7;
  typename Space::Elastic currLayer;
  typename Space::Elastic nextLayer;
  typename Space::Elastic tempLayer;
};

template<typename Space>
struct ElasticPointT<Space, 1>
{
  static const typename Space::IndexType paramCount = 7;
  typename Space::Elastic currLayer;
  typename Space::Elastic nextLayer;
  typename Space::Elastic subLayers[2];
};

template<typename Space>
struct ContactPointPairInfo
{
  typename Space::Vector3 normal;
  typename Space::IndexType typeIndex;
};

template<typename Space, int InterpolationOrder, typename InterpolatorT, int splitScheme>
class GcMeshT : public
  GeomMesh<typename Space::GeomSpace>
{
public:
  typedef Space                       ElasticSpace;
  typedef typename Space::GeomSpace   GeomSpace;
  typedef typename Space::Vector3     Vector3;
  typedef typename Space::Scalar      Scalar;
  typedef GcMeshT<Space, InterpolationOrder, InterpolatorT, splitScheme>
    GcMesh;

  typedef GeomMesh<GeomSpace> GeomMesh;

  const static int edgePointsCount = (InterpolationOrder - 1);
  const static int facePointsCount = (InterpolationOrder - 1) * (InterpolationOrder - 2) / 2;
  const static int cellPointsCount = (InterpolationOrder - 1) * (InterpolationOrder - 2) * (InterpolationOrder - 3) / 6;

  typedef ElasticPointT<Space, splitScheme>                             ElasticPoint;
  typedef typename Space::Elastic                                       Elastic;
  typedef typename Space::ElasticProperties                             ElasticProperties;

  struct ElasticPair
  {
    ElasticPoint *elasticPoints[2];
  };

  struct Vector3Pair
  {
    Vector3 vectors[2];
  };

  const static IndexType interpolationOrder = InterpolationOrder;

  std::vector<ElasticPoint> nodePoints;
  std::vector<ElasticPoint> edgePoints;
  std::vector<ElasticPoint> facePoints;
  std::vector<ElasticPoint> cellPoints;

  std::vector<char>         cellAxes;

  typedef InterpolatorT Interpolator;

  GcMeshT()
  {
  }

  ~GcMeshT()
  {
  }

  void LoadGeom(Vector3 *vertexPositions, IndexType *cellIndices, IndexType verticesCount, IndexType indicesCount,
                IndexType *submeshNodesCount, IndexType submeshesCount,
                FacePairIndices<GeomSpace> *contactFaces, IndexType *contactFacesCount, IndexType contactTypesCount,
                BoundaryFace<GeomSpace> *boundaryFaces, IndexType *boundaryFacesCount, IndexType boundaryTypesCount)
  {
    GeomMesh::LoadGeom(vertexPositions, cellIndices, verticesCount, indicesCount,
                       submeshNodesCount, submeshesCount,
                       contactFaces, contactFacesCount, contactTypesCount,
                       boundaryFaces, boundaryFacesCount, boundaryTypesCount);

    nodePoints.resize(this->GetNodesCount() * 1);
    edgePoints.resize(this->GetEdgesCount() * edgePointsCount);
    facePoints.resize(this->GetFacesCount() * facePointsCount);
    cellPoints.resize(this->GetCellsCount() * cellPointsCount);
    cellAxes  .resize(this->GetCellsCount());

    ComputeCellAxes();
  }

  IndexType GetTotalPointsCount()
  {
    IndexType totalPointsCount = 0;
    totalPointsCount += this->GetNodesCount() * 1;
    totalPointsCount += this->GetEdgesCount() * edgePointsCount;
    totalPointsCount += this->GetFacesCount() * facePointsCount;
    totalPointsCount += this->GetCellsCount() * cellPointsCount;
    return totalPointsCount;
  }

  IndexType GetTotalSubmeshPointsCount(IndexType submeshIndex)
  {
    IndexType totalPointsCount = 0;
    totalPointsCount += this->GetSubmeshNodesCount(submeshIndex) * 1;
    totalPointsCount += this->GetSubmeshEdgesCount(submeshIndex) * edgePointsCount;
    totalPointsCount += this->GetSubmeshFacesCount(submeshIndex) * facePointsCount;
    totalPointsCount += this->GetSubmeshCellsCount(submeshIndex) * cellPointsCount;
    return totalPointsCount;
  }


  void LoadState(IniStateMaker<Space> *stateMaker)
  {
    for(IndexType submeshIndex = 0; submeshIndex < this->submeshesCount; submeshIndex++)
    {
      IndexType elasticPointsCount = this->GetTotalSubmeshPointsCount(submeshIndex);

      UniversalGcIterator<GcMesh> elasticPointsStart(this, 0, submeshIndex);
      UniversalGcIterator<GcMesh> elasticPointsEnd(this, elasticPointsCount, submeshIndex);

      for(IndexType i = 0; i < elasticPointsCount; i++)
      {
        (elasticPointsStart + i).currValue() = stateMaker->GetValue((elasticPointsStart + i).currPos(), 1, 1, 1/*submeshIndex*/);
      }
    }
  }


  void MakeSnapshot3d(Elastic *destData, IndexType width, IndexType height, IndexType depth, Vector3 boxPoint1, Vector3 boxPoint2)
  {
    IndexType cells = this->GetCellsCount();

    Vector3 stepSize  (
                        (boxPoint2.x - boxPoint1.x) / (width - 1),
                        (boxPoint2.y - boxPoint1.y) / (height - 1),
                        (boxPoint2.z - boxPoint1.z) / (depth - 1)
                      );

    for(IndexType i = 0; i < cells; i++)
    {
      typename GeomMesh::Cell *cell = this->GetCell(i);
      typename GeomMesh::Vector3 points[4];
      points[0] = this->GetNode(cell->incidentNodes[0])->pos;
      points[1] = this->GetNode(cell->incidentNodes[1])->pos;
      points[2] = this->GetNode(cell->incidentNodes[2])->pos;
      points[3] = this->GetNode(cell->incidentNodes[3])->pos;
      Scalar 
        xmin = points[0].x, xmax = points[0].x,
        ymin = points[0].y, ymax = points[0].y,
        zmin = points[0].z, zmax = points[0].z;
      for(int j = 0; j < 4; j++)
      {
        if(points[j].x > 99.0 && points[j].y < 1.0 && points[j].z < 1.0)
        {
          int pp = 1;
        }
        if(points[j].x < xmin) xmin = points[j].x;
        if(points[j].x > xmax) xmax = points[j].x;
        if(points[j].y < ymin) ymin = points[j].y;
        if(points[j].y > ymax) ymax = points[j].y;
        if(points[j].z < zmin) zmin = points[j].z;
        if(points[j].z > zmax) zmax = points[j].z;
      }
      IndexType ixmin = IndexType(((xmin - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * Scalar(width  - 1));
      IndexType ixmax = IndexType(((xmax - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * Scalar(width  - 1)) + 1;
      IndexType iymin = IndexType(((ymin - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * Scalar(height - 1));
      IndexType iymax = IndexType(((ymax - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * Scalar(height - 1)) + 1;
      IndexType izmin = IndexType(((zmin - boxPoint1.z) / (boxPoint2.z - boxPoint1.z)) * Scalar(depth  - 1));
      IndexType izmax = IndexType(((zmax - boxPoint1.z) / (boxPoint2.z - boxPoint1.z)) * Scalar(depth  - 1)) + 1;

      if(ixmin < 0) ixmin = 0;
      if(ixmin > width - 1) ixmin = width - 1;
      if(ixmax < 0) ixmax = 0;
      if(ixmax > width - 1) ixmax = width - 1;

      if(iymin < 0) iymin = 0;
      if(iymin > height - 1) iymin = height - 1;
      if(iymax < 0) iymax = 0;
      if(iymax > height - 1) iymax = height - 1;

      if(izmin < 0) izmin = 0;
      if(izmin > depth - 1) izmin = depth - 1;
      if(izmax < 0) izmax = 0;
      if(izmax > depth - 1) izmax = depth - 1;

      for(IndexType x = ixmin; x <= ixmax; x++)
      {
        for(IndexType y = iymin; y <= iymax; y++)
        {
          for(IndexType z = izmin; z <= izmax; z++)
          {
            typename GeomMesh::Vector3 point(boxPoint1.x + (boxPoint2.x - boxPoint1.x) * (Scalar(x) / Scalar(width - 1)),
                                boxPoint1.y + (boxPoint2.y - boxPoint1.y) * (Scalar(y) / Scalar(height - 1)),
                                boxPoint1.z + (boxPoint2.z - boxPoint1.z) * (Scalar(z) / Scalar(depth - 1)));
            if(PointInCell<GeomSpace>(points, point))
            {
              destData[z * width + y * width * depth + x] = interpolateElastic(i, point);
//              destData[z * width + y * width * depth + x].v = Vector3(1, 0, 0);
            }
          }
        }
      }
    }
  }

  void MakeBoundarySnapshot3d(Elastic *destData, IndexType width, IndexType height, IndexType depth, Vector3 boxPoint1, Vector3 boxPoint2)
  {
    IndexType cells = this->GetCellsCount();

    ContactGroupGcIterator<ElasticMesh> contactPointsStart(this, 0);
    ContactGroupGcIterator<ElasticMesh> contactPointsEnd = contactPointsStart + GetContactGroupPointsCount();

    /*for(IndexType x = 0; x < width; x++)
    {
      for(IndexType y = 0; y < height; y++)
      {
        for(IndexType z = 0; z < depth; z++)
        {
          destData[z * width + y * width * depth + x].v = Vector3(0, 0, 0);
        }
      }
    }*/

    for(IndexType pointIndex = 0; pointIndex < GetContactGroupPointsCount(); pointIndex++)
    {
      IndexType sidesCount = (contactPointsStart + pointIndex).GetContactSidesCount();
      if(sidesCount == 2)
      {
        Vector3 pos = ((contactPointsStart + pointIndex).currPos(0) + (contactPointsStart + pointIndex).currPos(1)) * Scalar(0.5);
        //Vector3 pos = (contactPointsStart + pointIndex).currPos(0);
        IndexType type = (contactPointsStart + pointIndex).GetContactTypeIndex(0, 1);
        IndexType x = IndexType(((pos.x - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * width + Scalar(0.5));
        IndexType y = IndexType(((pos.y - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * height + Scalar(0.5));
        IndexType z = IndexType(((pos.z - boxPoint1.z) / (boxPoint2.z - boxPoint1.z)) * depth + Scalar(0.5));

        if(x > 0 && x < width && y > 0 && y < height && z > 0 && z < depth)
        {
          switch(type)
          {
            case IndexType(-1):
            {
              destData[z * width + y * width * depth + x].v.x = 0.5;
            }
            break;
            case 0:
            {
              destData[z * width + y * width * depth + x].v.x = 0.1;
            }
            break;
            case 1:
            {
              destData[z * width + y * width * depth + x].v.x = 0.1;
            }
            break;
            case 2:
            {
              destData[z * width + y * width * depth + x].v.x = 0.1;
            }
            break;
          }
        }
      }

      if(sidesCount == 1)
      {
        Vector3 pos = (contactPointsStart + pointIndex).currPos(0);
        IndexType type = (contactPointsStart + pointIndex).GetContactTypeIndex(0, -1);
        IndexType x = IndexType(((pos.x - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * width + Scalar(0.5));
        IndexType y = IndexType(((pos.y - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * height + Scalar(0.5));
        IndexType z = IndexType(((pos.z - boxPoint1.z) / (boxPoint2.z - boxPoint1.z)) * depth + Scalar(0.5));

        if(x > 0 && x < width && y > 0 && y < height && z > 0 && z < depth)
        {
          destData[z * width + y * width * depth + x].v = (contactPointsStart + pointIndex).currNormal(0, -1);
          /*switch(type)
          {
            case IndexType(-1):
            {
              destData[z * width + y * width * depth + x].v.x = 0.5;
            }
            break;
            case 0:
            {
              destData[z * width + y * width * depth + x].v.x = 0.3;
            }
            break;
            case 1:
            {
              destData[z * width + y * width * depth + x].v.x = 0.2;
            }
            break;
            case 2:
            {
              destData[z * width + y * width * depth + x].v.x = 0.1;
            }
            break;
          }*/
        }
      }
    }
  }

  template<IndexType paramIndex>
  Scalar interpolateInCell(IndexType cellIndex, Vector3 point) const
  {

    const IndexType dataSize = (InterpolationOrder + 1) * (InterpolationOrder + 2) * (InterpolationOrder + 3) / 6;
    Scalar interpolatorData[dataSize];
    Vector3 points[5];

    points[0] = point;

    IndexType offset = 0;

    //nodes
    //IndexType nodeOffset = offset;
    for(IndexType i = 0; i < 4; i++)
    {
//      Scalar p = nodeIt.GetNode()->data.currLayer(paramIndex);
      IndexType nodeIndex = cells[cellIndex].incidentNodes[i];
      Scalar p = GetNodeElasticData(nodeIndex)->currLayer(paramIndex);
      points[i + 1] = nodes[nodeIndex].pos;
      interpolatorData[offset++] = p;
    }
    
    {
      typename GeomMesh::Vector3 points[4];
      points[0] = this->GetNode(this->GetCell(cellIndex)->incidentNodes[0])->pos;
      points[1] = this->GetNode(this->GetCell(cellIndex)->incidentNodes[1])->pos;
      points[2] = this->GetNode(this->GetCell(cellIndex)->incidentNodes[2])->pos;
      points[3] = this->GetNode(this->GetCell(cellIndex)->incidentNodes[3])->pos;

      if(!PointInCell<GeomSpace>(points, point))
      {
        printf("Point in cell interpolation does not belong to the cell\n ");
        /*printf("Point0 : %f %f %f\n", points[0].x, points[0].y, points[0].z);
        printf("Point1 : %f %f %f\n", points[1].x, points[1].y, points[1].z);
        printf("Point2 : %f %f %f\n", points[2].x, points[2].y, points[2].z);
        printf("Point3 : %f %f %f\n", points[3].x, points[3].y, points[3].z);
        printf("Point  : %f %f %f\n", point.x    , point.y    , point.z    );
        while(1);*/
      }
    }

    if(edgePointsCount > 0)
    {
      //IndexType edgeOffset = offset;
      for(IndexType i = 0; i < 6; i++)
      {
        IndexType edgeIndex = cells[cellIndex].incidentEdges[i];
        IndexType controlPointsCount = edgePointsCount;
        for(IndexType j = 0; j < controlPointsCount; j++)
        {
          Scalar p = GetEdgeElasticData(edgeIndex, j)->currLayer(paramIndex);
          interpolatorData[offset++] = p;
        }
      }
    }

    if(facePointsCount > 0)
    {
      //IndexType faceOffset = offset;
      for(IndexType i = 0; i < 4; i++)
      {
        IndexType faceIndex = cells[cellIndex].incidentFaces[i];
        IndexType controlPointsCount = facePointsCount;
        for(IndexType j = 0; j < controlPointsCount; j++)
        {
          Scalar p = GetFaceElasticData(faceIndex, j)->currLayer(paramIndex);
          interpolatorData[offset++] = p;
        }
      }
    }

    if(cellPointsCount > 0)
    {
      //IndexType cellOffset = offset;
      {
        IndexType controlPointsCount = cellPointsCount;
        for(IndexType j = 0; j < controlPointsCount; j++)
        {
          Scalar p = this->GetCellElasticData(cellIndex, j)->currLayer(paramIndex);
          interpolatorData[offset++] = p;
        }
      }
    }
     
    if(offset != dataSize)
    {
      //int gg = 1;
      //something's wrong
    }

    return Interpolator::Interpolate( points[0], 
                                      points + 1,
                                      interpolatorData, cellAxes[cellIndex]);
  }
  Elastic interpolateElastic(IndexType cellIndex, Vector3 point) const
  {
    Elastic m;
    //velocity
    m(0) = interpolateInCell<0>(cellIndex, point);
    m(1) = interpolateInCell<1>(cellIndex, point);
    m(2) = interpolateInCell<2>(cellIndex, point);

    //sigma
    m(3) = interpolateInCell<3>(cellIndex, point);
    m(4) = interpolateInCell<4>(cellIndex, point);
    m(5) = interpolateInCell<5>(cellIndex, point);
    m(6) = interpolateInCell<6>(cellIndex, point);
    m(7) = interpolateInCell<7>(cellIndex, point);
    m(8) = interpolateInCell<8>(cellIndex, point);
    return m;
  }

  static void RotateEdgeData(Elastic *edgePoints, IndexType *targetIndices)
  {
  }

  static void RotateFaceData(Elastic *facePoints, IndexType *targetIndices)
  {
  }

  static void RotateCellData(Elastic *cellPoints, IndexType *targetIndices)
  {
  }
/*  static ElasticFaceData RotateFaceData(const ElasticFaceData &srcMeshData, IndexType *targetIndices)
  {
    ElasticFaceData dataToRotate = srcMeshData;
    //Interpolator::R
    return dataToRotate;
  }

  static ElasticCellData RotateCellData(const ElasticCellData &srcMeshData, IndexType *targetIndices)
  {
    ElasticCellData dataToRotate = srcMeshData;
    //Interpolator::R
    return dataToRotate;
  }*/

  IndexType GetSubmeshNodeGlobalIndex(IndexType submeshNodeIndex, IndexType submeshIndex)
  {
    return submeshNodeIndex + this->submeshInfos[submeshIndex].firstNodeIndex;
  }
  IndexType GetSubmeshEdgeGlobalIndex(IndexType submeshEdgeIndex, IndexType submeshIndex)
  {
    return submeshEdgeIndex + this->submeshInfos[submeshIndex].firstEdgeIndex;
  }
  IndexType GetSubmeshFaceGlobalIndex(IndexType submeshFaceIndex, IndexType submeshIndex)
  {
    return submeshFaceIndex + this->submeshInfos[submeshIndex].firstFaceIndex;
  }
  IndexType GetSubmeshCellGlobalIndex(IndexType submeshCellIndex, IndexType submeshIndex)
  {
    return submeshCellIndex + this->submeshInfos[submeshIndex].firstCellIndex;
  }

  IndexType GetSubmeshNodePointGlobalIndex(IndexType submeshNodePointIndex, IndexType submeshIndex)
  {
    return submeshNodePointIndex + this->submeshInfos[submeshIndex].firstNodeIndex;
  }
  IndexType GetSubmeshEdgePointGlobalIndex(IndexType submeshEdgePointIndex, IndexType submeshIndex)
  {
    return submeshEdgePointIndex + this->submeshInfos[submeshIndex].firstEdgeIndex * edgePointsCount;
  }
  IndexType GetSubmeshFacePointGlobalIndex(IndexType submeshFacePointIndex, IndexType submeshIndex)
  {
    return submeshFacePointIndex + this->submeshInfos[submeshIndex].firstFaceIndex * facePointsCount;
  }
  IndexType GetSubmeshCellPointGlobalIndex(IndexType submeshCellPointIndex, IndexType submeshIndex)
  {
    return submeshCellPointIndex + this->submeshInfos[submeshIndex].firstCellIndex * cellPointsCount;
  }

  Vector3 GetNodePointPos(IndexType nodeIndex)
  {
    return this->GetNode(nodeIndex)->pos;
  }

  Vector3 GetSubmeshNodePointPos(IndexType submeshNodeIndex, IndexType submeshIndex)
  {
    return this->GetNodePointPos(this->GetSubmeshNodeGlobalIndex(submeshNodeIndex, submeshIndex));
  }

  Vector3 GetEdgePointPos(IndexType edgeIndex, IndexType pointNum)
  {
    return Interpolator::GetEdgePointPos( this->GetNode(this->GetEdge(edgeIndex)->incidentNodes[0])->pos,
                                          this->GetNode(this->GetEdge(edgeIndex)->incidentNodes[1])->pos, pointNum);
  }
  Vector3 GetSubmeshEdgePointPos(IndexType submeshEdgeIndex, IndexType pointNum, IndexType submeshIndex)
  {
    return this->GetEdgePointPos(this->GetSubmeshEdgeGlobalIndex(submeshEdgeIndex, submeshIndex), pointNum);
  }

  Vector3 GetEdgePointPos(IndexType globalEdgePointIndex)
  {
    if(edgePointsCount > 0) 
      return GetEdgePointPos(globalEdgePointIndex / edgePointsCount, globalEdgePointIndex % edgePointsCount);
    else
      return Vector3(0, 0, 0);
  }
  Vector3 GetSubmeshEdgePointPos(IndexType submeshEdgePointIndex, IndexType submeshIndex)
  {
    return this->GetEdgePointPos(this->GetSubmeshEdgePointGlobalIndex(submeshEdgePointIndex, submeshIndex));
  }

  Vector3 GetFacePointPos(IndexType faceIndex, IndexType pointNum)
  {
    return Interpolator::GetFacePointPos( this->GetNode(this->GetFace(faceIndex)->incidentNodes[0])->pos,
                                          this->GetNode(this->GetFace(faceIndex)->incidentNodes[1])->pos,
                                          this->GetNode(this->GetFace(faceIndex)->incidentNodes[2])->pos, pointNum);
  }
  Vector3 GetSubmeshFacePointPos(IndexType submeshFaceIndex, IndexType pointNum, IndexType submeshIndex)
  {
    return this->GetFacePointPos(this->GetSubmeshFaceGlobalIndex(submeshFaceIndex, submeshIndex), pointNum);
  }

  Vector3 GetFacePointPos(IndexType globalFacePointIndex)
  {
    if(facePointsCount > 0)
      return GetFacePointPos(globalFacePointIndex / facePointsCount, globalFacePointIndex % facePointsCount);
    else
      return Vector3(0, 0, 0);
  }
  Vector3 GetSubmeshFacePointPos(IndexType submeshFacePointIndex, IndexType submeshIndex)
  {
    return this->GetFacePointPos(this->GetSubmeshFacePointGlobalIndex(submeshFacePointIndex, submeshIndex));
  }

  Vector3 GetCellPointPos(IndexType cellIndex, IndexType pointNum)
  {
    return Interpolator::GetCellPointPos( this->GetNode(this->GetCell(cellIndex)->incidentNodes[0])->pos,
                                          this->GetNode(this->GetCell(cellIndex)->incidentNodes[1])->pos,
                                          this->GetNode(this->GetCell(cellIndex)->incidentNodes[2])->pos,
                                          this->GetNode(this->GetCell(cellIndex)->incidentNodes[3])->pos, pointNum);
  }
  Vector3 GetSubmeshCellPointPos(IndexType submeshCellIndex, IndexType pointNum, IndexType submeshIndex)
  {
    return this->GetCellPointPos(this->GetSubmeshCellGlobalIndex(submeshCellIndex, submeshIndex), pointNum);
  }

  Vector3 GetCellPointPos(IndexType globalCellPointIndex)
  {
    if(cellPointsCount > 0)
      return GetCellPointPos(globalCellPointIndex / cellPointsCount, globalCellPointIndex % cellPointsCount);
    else
      return Vector3(0, 0, 0);
  }
  Vector3 GetSubmeshCellPointPos(IndexType submeshCellPointIndex, IndexType submeshIndex)
  {
    return this->GetCellPointPos(this->GetSubmeshCellPointGlobalIndex(submeshCellPointIndex, submeshIndex));
  }

  ElasticPoint *GetNodeElasticData(IndexType nodeIndex)
  {
    return &(nodePoints[nodeIndex]);
  }
  ElasticPoint *GetSubmeshNodeElasticData(IndexType submeshNodeIndex, IndexType submeshIndex)
  {
    return this->GetNodeElasticData(this->GetSubmeshNodePointGlobalIndex(submeshNodeIndex, submeshIndex));
  }

  ElasticPoint *GetEdgeElasticData(IndexType edgeIndex, IndexType pointIndex)
  {
    return GetEdgeElasticData(edgeIndex * edgePointsCount + pointIndex);
  }
  ElasticPoint *GetFaceElasticData(IndexType faceIndex, IndexType pointIndex)
  {
    return GetFaceElasticData(faceIndex * facePointsCount + pointIndex);
  }
  ElasticPoint *GetCellElasticData(IndexType cellIndex, IndexType pointIndex)
  {
    return GetCellElasticData(cellIndex * cellPointsCount + pointIndex);
  }

  ElasticPoint *GetEdgeElasticData(IndexType globalPointIndex)
  {
    return &(edgePoints[globalPointIndex]);
  }
  ElasticPoint *GetSubmeshEdgeElasticData(IndexType submeshEdgePointIndex, IndexType submeshIndex)
  {
    return this->GetEdgeElasticData(this->GetSubmeshEdgePointGlobalIndex(submeshEdgePointIndex, submeshIndex));
  }

  ElasticPoint *GetFaceElasticData(IndexType globalPointIndex)
  {
    return &(facePoints[globalPointIndex]);
  }
  ElasticPoint *GetSubmeshFaceElasticData(IndexType submeshFacePointIndex, IndexType submeshIndex)
  {
    return this->GetFaceElasticData(this->GetSubmeshFacePointGlobalIndex(submeshFacePointIndex, submeshIndex));
  }

  ElasticPoint *GetCellElasticData(IndexType globalPointIndex)
  {
    return &(cellPoints[globalPointIndex]);
  }
  ElasticPoint *GetSubmeshCellElasticData(IndexType submeshCellPointIndex, IndexType submeshIndex)
  {
    return this->GetCellElasticData(this->GetSubmeshCellPointGlobalIndex(submeshCellPointIndex, submeshIndex));
  }

  const ElasticPoint *GetNodeElasticData(IndexType nodeIndex) const
  {
    return &(nodePoints[nodeIndex]);
  }

  const ElasticPoint *GetEdgeElasticData(IndexType edgeIndex, IndexType pointIndex) const
  {
    return GetEdgeElasticData(edgeIndex * edgePointsCount + pointIndex);
  }
  const ElasticPoint *GetFaceElasticData(IndexType faceIndex, IndexType pointIndex) const
  {
    return GetFaceElasticData(faceIndex * facePointsCount + pointIndex);
  }
  const ElasticPoint *GetCellElasticData(IndexType cellIndex, IndexType pointIndex) const
  {
    return GetCellElasticData(cellIndex * cellPointsCount + pointIndex);
  }

  const ElasticPoint *GetEdgeElasticData(IndexType globalPointIndex) const
  {
    return &(edgePoints[globalPointIndex]);
  }
  const ElasticPoint *GetFaceElasticData(IndexType globalPointIndex) const
  {
    return &(facePoints[globalPointIndex]);
  }
  const ElasticPoint *GetCellElasticData(IndexType globalPointIndex) const
  {
    return &(cellPoints[globalPointIndex]);
  }


  Vector3 GetBorderFacePointNormal(IndexType faceBorderIndex, IndexType pointNum)
  {
    IndexType faceIndex = this->GetBorderFaceIndex(faceBorderIndex);

    IndexType faceNode0 = this->faces[faceIndex].incidentNodes[0];
    IndexType faceNode1 = this->faces[faceIndex].incidentNodes[1];
    IndexType faceNode2 = this->faces[faceIndex].incidentNodes[2];

    IndexType faceBorderNode0 = this->nodeGlobalToBorder[faceNode0];
    IndexType faceBorderNode1 = this->nodeGlobalToBorder[faceNode1];
    IndexType faceBorderNode2 = this->nodeGlobalToBorder[faceNode2];

    Vector3 nodeNormal0 = this->GetBorderNodeNormal(faceBorderNode0);
    Vector3 nodeNormal1 = this->GetBorderNodeNormal(faceBorderNode1);
    Vector3 nodeNormal2 = this->GetBorderNodeNormal(faceBorderNode2);

    Vector3 normal = Interpolator::GetFacePointPos( nodeNormal0, nodeNormal1, nodeNormal2, pointNum);
    return normal * (Scalar(1.0) / absolute(normal));
  }

  Vector3 GetContactEdgePointNormal(IndexType edgeContactIndex, IndexType pointNum)
  {
    IndexType edgeIndex = this->GetContactEdgePair(edgeContactIndex).edges[0];

    IndexType edgeNode0 = this->edges[edgeIndex].incidentNodes[0];
    IndexType edgeNode1 = this->edges[edgeIndex].incidentNodes[1];

    IndexType edgeContactNode0 = this->nodeGlobalToContact[edgeNode0];
    IndexType edgeContactNode1 = this->nodeGlobalToContact[edgeNode1];

    Vector3 nodeNormal0 = this->GetContactNodeNormals(edgeContactNode0).normals[0];
    Vector3 nodeNormal1 = this->GetContactNodeNormals(edgeContactNode1).normals[0];

    Vector3 normal = Interpolator::GetEdgePointPos( nodeNormal0, nodeNormal1, pointNum);
    return normal * (Scalar(1.0) / absolute(normal));
  }

  Vector3 GetContactFacePointNormal(IndexType faceContactIndex, IndexType pointNum)
  {
    IndexType faceIndex = this->GetContactFacePair(faceContactIndex).faces[0];

    IndexType faceNode0 = this->faces[faceIndex].incidentNodes[0];
    IndexType faceNode1 = this->faces[faceIndex].incidentNodes[1];
    IndexType faceNode2 = this->faces[faceIndex].incidentNodes[2];

    IndexType faceContactNode0 = this->nodeGlobalToContact[faceNode0];
    IndexType faceContactNode1 = this->nodeGlobalToContact[faceNode1];
    IndexType faceContactNode2 = this->nodeGlobalToContact[faceNode2];

    Vector3 nodeNormal0 = this->GetContactNodeNormals(faceContactNode0).normals[0];
    Vector3 nodeNormal1 = this->GetContactNodeNormals(faceContactNode1).normals[0];
    Vector3 nodeNormal2 = this->GetContactNodeNormals(faceContactNode2).normals[0];

    Vector3 normal = Interpolator::GetFacePointPos( nodeNormal0, nodeNormal1, nodeNormal2, pointNum);
    return normal * (Scalar(1.0) / absolute(normal));
  }

  Vector3 GetContactNodeGroupPointPos(IndexType groupIndex, IndexType sideIndex)
  {
    ContactNode<GeomSpace> contactNode = this->GetContactNodeInGroup(groupIndex, sideIndex);
    return GetNode(contactNode.nodeIndex)->pos;
  }

  Vector3 GetContactEdgeGroupPointPos(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex)
  {
    ContactEdge<GeomSpace> referenceContactEdge = this->GetContactEdgeInGroup(groupIndex, 0);
    ContactEdge<GeomSpace> contactEdge = this->GetContactEdgeInGroup(groupIndex, sideIndex);
    IndexType orientation = GetRelativeOrientationEdge<GeomSpace>(referenceContactEdge.incidentNodes, contactEdge.incidentNodes);

    IndexType transposedPointIndex = GetCorrespondingEdgeIndex(orientation, pointIndex, InterpolationOrder);

    return GetEdgePointPos(contactEdge.edgeIndex, transposedPointIndex);
  }

  Vector3 GetContactFaceGroupPointPos(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex)
  {
    ContactFace<GeomSpace> referenceContactFace = this->GetContactFaceInGroup(groupIndex, 0);
    ContactFace<GeomSpace> contactFace = this->GetContactFaceInGroup(groupIndex, sideIndex);
    IndexType orientation = GetRelativeOrientationSide<GeomSpace>(referenceContactFace.incidentNodes, contactFace.incidentNodes);

    IndexType transposedPointIndex = GetCorrespondingFaceIndex(orientation, pointIndex, InterpolationOrder);

    return GetFacePointPos(contactFace.faceIndex, transposedPointIndex);
  }

  ElasticPoint *GetContactNodeGroupElasticData(IndexType groupIndex, IndexType sideIndex)
  {
    ContactNode<GeomSpace> contactNode = this->GetContactNodeInGroup(groupIndex, sideIndex);
    return GetNodeElasticData(contactNode.nodeIndex);
  }

  ElasticPoint *GetContactEdgeGroupElasticData(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex)
  {
    ContactEdge<GeomSpace> referenceContactEdge = this->GetContactEdgeInGroup(groupIndex, 0);
    ContactEdge<GeomSpace> contactEdge = this->GetContactEdgeInGroup(groupIndex, sideIndex);
    IndexType orientation = GetRelativeOrientationEdge<GeomSpace>(referenceContactEdge.incidentNodes, contactEdge.incidentNodes);
    IndexType transposedPointIndex = GetCorrespondingEdgeIndex(orientation, pointIndex, InterpolationOrder);

    return GetEdgeElasticData(contactEdge.edgeIndex, transposedPointIndex);
  }

  ElasticPoint *GetContactFaceGroupElasticData(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex)
  {
    ContactFace<GeomSpace> referenceContactFace = this->GetContactFaceInGroup(groupIndex, 0);
    ContactFace<GeomSpace> contactFace = this->GetContactFaceInGroup(groupIndex, sideIndex);
    IndexType orientation = GetRelativeOrientationSide<GeomSpace>(referenceContactFace.incidentNodes, contactFace.incidentNodes);
    IndexType transposedPointIndex = GetCorrespondingFaceIndex(orientation, pointIndex, InterpolationOrder);

    return GetFaceElasticData(contactFace.faceIndex, transposedPointIndex);
  }

  /*Vector3 GetContactNodeGroupNormal(IndexType groupIndex, IndexType sideIndex0, IndexType sideIndex1)
  {
    ContactInfo<Vector3> contactInfo = this->GetContactNodeGroupInfo(groupIndex, sideIndex0, sideIndex1);
    return contactInfo.normal;
  }*/

  /*ContactInfo<Vector3> GetNodesContactInfo(IndexType nodeIndex0, IndexType nodeIndex1)
  {
    IndexType groupIndex = this->nodes[nodeIndex0].contactGroupIndex;
    assert(groupIndex == this->nodes[nodeIndex1].contactGroupIndex);

    IndexType sideIndex0 = this->GetContactSideInGroup(groupIndex, nodeIndex0);
    IndexType sideIndex1 = this->GetContactSideInGroup(groupIndex, nodeIndex1);

    assert(sideIndex0 != -1);
    assert(sideIndex1 != -1);

    return this->GetContactNodeGroupInfo(groupIndex, sideIndex0, sideIndex1);
  }

  ContactInfo<Vector3> GetNodeBorderInfo(IndexType nodeIndex)
  {
    IndexType groupIndex = this->nodes[nodeIndex].contactGroupIndex;
    IndexType sideIndex = this->GetContactSideInGroup(groupIndex, nodeIndex);
    assert(sideIndex != -1);
    return this->GetContactNodeGroupInfo(groupIndex, sideIndex, -1);
  }

  ContactInfo<Vector3> GetNodeFreeBorderInfo(IndexType nodeIndex)
  {
    IndexType groupIndex = this->nodes[nodeIndex].contactGroupIndex;
    return this->GetContactNodeGroupInfo(groupIndex, -1, -1);
  }

  ContactInfo<Vector3> GetContactEdgeGroupInfo(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex0, IndexType sideIndex1)
  {
    ContactInfo<Vector3> infos[2];
    ContactInfo<Vector3> resultInfo;

    ContactEdge<Vector3> contactEdge0;
    ContactEdge<Vector3> contactEdge1;
    if((sideIndex0 != -1) && (sideIndex1 != -1))
    {
      contactEdge0 = this->GetContactEdgeInGroup(groupIndex, sideIndex0);
      contactEdge1 = this->GetContactEdgeInGroup(groupIndex, sideIndex1);
      infos[0] = GetNodesContactInfo(contactEdge0.incidentNodes[0], contactEdge1.incidentNodes[0]);
      infos[1] = GetNodesContactInfo(contactEdge0.incidentNodes[1], contactEdge1.incidentNodes[1]);
    }else
    if((sideIndex0 != -1) && (sideIndex1 == -1))
    {
      contactEdge0 = this->GetContactEdgeInGroup(groupIndex, sideIndex0);
      infos[0] = GetNodeBorderInfo(contactEdge0.incidentNodes[0]);
      infos[1] = GetNodeBorderInfo(contactEdge0.incidentNodes[1]);
    }else
    if((sideIndex0 == -1) && (sideIndex1 != -1))
    {
      contactEdge1 = this->GetContactEdgeInGroup(groupIndex, sideIndex1);
      infos[0] = GetNodeBorderInfo(contactEdge1.incidentNodes[0]);
      infos[1] = GetNodeBorderInfo(contactEdge1.incidentNodes[1]);
    }else
    {
      contactEdge0 = this->GetContactEdgeInGroup(groupIndex, 0);
      infos[0] = GetNodeFreeBorderInfo(contactEdge0.incidentNodes[0]);
      infos[1] = GetNodeFreeBorderInfo(contactEdge0.incidentNodes[1]);
    }

    Vector3 normal = Vector3(0, 0, 0);
    resultInfo.typeIndex = IndexType(-1);
    for(IndexType nodeNumber = 0; nodeNumber < 2; nodeNumber++)
    {
      if(infos[nodeNumber].typeIndex != IndexType(-1))
      {
         if((infos[nodeNumber].typeIndex > resultInfo.typeIndex) || (resultInfo.typeIndex == IndexType(-1)))
         {
           resultInfo.typeIndex = infos[nodeNumber].typeIndex;
           normal = Vector3(0, 0, 0);
         }
         if(infos[nodeNumber].typeIndex == resultInfo.typeIndex)
         {
           normal += infos[nodeNumber].normal;
         }
      }
    }

    Scalar len = absolute(normal);
    if(len > 0.0)
      normal = normal / len;
    else
    {
      normal = Vector3(Scalar(1.0), 0, 0);
    }
    resultInfo.normal = normal;
    return resultInfo;
  }

  ContactInfo<Vector3> GetContactFaceGroupInfo(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex0, IndexType sideIndex1)
  {
    ContactFace<Vector3> contactFace0;
    ContactFace<Vector3> contactFace1;

    ContactInfo<Vector3> infos[3];
    ContactInfo<Vector3> resultInfo;

    if((sideIndex0 != -1) && (sideIndex1 != -1))
    {
      contactFace0 = this->GetContactFaceInGroup(groupIndex, sideIndex0);
      contactFace1 = this->GetContactFaceInGroup(groupIndex, sideIndex1);

      infos[0] = GetNodesContactInfo(contactFace0.incidentNodes[0], contactFace1.incidentNodes[0]);
      infos[1] = GetNodesContactInfo(contactFace0.incidentNodes[1], contactFace1.incidentNodes[1]);
      infos[2] = GetNodesContactInfo(contactFace0.incidentNodes[2], contactFace1.incidentNodes[2]);
    }else
    if((sideIndex0 != -1) && (sideIndex1 == -1))
    {
      contactFace0 = this->GetContactFaceInGroup(groupIndex, sideIndex0);
      infos[0] = GetNodeBorderInfo(contactFace0.incidentNodes[0]);
      infos[1] = GetNodeBorderInfo(contactFace0.incidentNodes[1]);
      infos[2] = GetNodeBorderInfo(contactFace0.incidentNodes[2]);
    }else
    if((sideIndex0 == -1) && (sideIndex1 != -1))
    {
      contactFace0 = this->GetContactFaceInGroup(groupIndex, sideIndex1);
      infos[0] = GetNodeBorderInfo(contactFace0.incidentNodes[0]);
      infos[1] = GetNodeBorderInfo(contactFace0.incidentNodes[1]);
      infos[2] = GetNodeBorderInfo(contactFace0.incidentNodes[2]);
    }else
    {
      contactFace0 = this->GetContactFaceInGroup(groupIndex, 0);
      infos[0] = GetNodeFreeBorderInfo(contactFace0.incidentNodes[0]);
      infos[1] = GetNodeFreeBorderInfo(contactFace0.incidentNodes[1]);
      infos[2] = GetNodeFreeBorderInfo(contactFace0.incidentNodes[2]);
    }



    Vector3 normal = Vector3(-1, 0, 0);
    resultInfo.typeIndex = IndexType(-1);
    for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
    {
      if(infos[nodeNumber].typeIndex != IndexType(-1))
      {
         if((infos[nodeNumber].typeIndex > resultInfo.typeIndex) || (resultInfo.typeIndex == IndexType(-1)))
         {
           resultInfo.typeIndex = infos[nodeNumber].typeIndex;
           normal = Vector3(0, 0, 0);
         }
         if(infos[nodeNumber].typeIndex == resultInfo.typeIndex)
         {
           normal += infos[nodeNumber].normal;
         }
      }
    }

    Scalar len = absolute(normal);
    if(len > 0.0)
      normal = normal / len;
    else
    {
      normal = Vector3(Scalar(1.0), 0, 0);
    }
    resultInfo.normal = normal;
    return resultInfo;
  }*/



  /*ContactNodePairInfo<Vector3> GetContactNodePairInfo(IndexType nodeIndex0, IndexType nodeIndex1)
  {
    IndexType groupIndex = this->nodes[nodeIndex0].contactGroupIndex;
    assert(groupIndex == this->nodes[nodeIndex1].contactGroupIndex);

    IndexType sideIndex0 = this->GetContactSideInGroup(groupIndex, nodeIndex0);
    IndexType sideIndex1 = this->GetContactSideInGroup(groupIndex, nodeIndex1);

    assert(sideIndex0 != -1);
    assert(sideIndex1 != -1);

    return this->GetContactNodePairInfo(groupIndex, sideIndex0, sideIndex1);
  }*/

  ContactPointPairInfo<GeomSpace> GetContactEdgePointPairInfo(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex0, IndexType sideIndex1)
  {
    ContactPointPairInfo<GeomSpace> resultInfo;
    ContactEdgePairInfo<GeomSpace> edgePairInfo = this->GetContactEdgePairInfo(groupIndex, sideIndex0, sideIndex1);

    resultInfo.typeIndex = edgePairInfo.typeIndex;
    resultInfo.normal = Vector3(-1, 0, 0);

    if(resultInfo.typeIndex != IndexType(-1))
    {
      Vector3 normal = Interpolator::GetEdgePointPos(edgePairInfo.normals[0], edgePairInfo.normals[1], pointIndex);
      resultInfo.normal = normal.GetNorm();
    }

    return resultInfo;
  }

  ContactPointPairInfo<GeomSpace> GetContactFacePointPairInfo(IndexType groupIndex, IndexType pointIndex, IndexType sideIndex0, IndexType sideIndex1)
  {
    ContactPointPairInfo<GeomSpace> resultInfo;

    ContactFacePairInfo<GeomSpace> facePairInfo = this->GetContactFacePairInfo(groupIndex, sideIndex0, sideIndex1);

    resultInfo.typeIndex = facePairInfo.typeIndex;
    resultInfo.normal = Vector3(-1, 0, 0);

    if(resultInfo.typeIndex != IndexType(-1))
    {
      Vector3 normal = Interpolator::GetFacePointPos( facePairInfo.normals[0], facePairInfo.normals[1], facePairInfo.normals[2], pointIndex);
      resultInfo.normal = normal.GetNorm();
    }

    return resultInfo;
  }

  
  IndexType GetContactGroupPointsCount()
  {
    IndexType totalPointsCount = 0;
    totalPointsCount += this->GetContactNodeGroupsCount() * 1;
    totalPointsCount += this->GetContactEdgeGroupsCount() * edgePointsCount;
    totalPointsCount += this->GetContactFaceGroupsCount() * facePointsCount;
    return totalPointsCount;
  }


  Scalar GetMinCellHeight()
  {
    Scalar minHeight = Scalar(1e5);
    IndexType cells = this->GetCellsCount();
    for(IndexType i = 0; i < cells; i++)
    {
      typename GeomMesh::Cell *cell = this->GetCell(i);
      typename GeomMesh::Vector3 points[4];
      points[0] = this->GetNode(cell->incidentNodes[0])->pos;
      points[1] = this->GetNode(cell->incidentNodes[1])->pos;
      points[2] = this->GetNode(cell->incidentNodes[2])->pos;
      points[3] = this->GetNode(cell->incidentNodes[3])->pos;

      Vector3 n0 = outer_product(points[1] - points[2], points[3] - points[2]);
      Scalar h0 = fabs(scalar_product(points[0] - points[2], n0) / absolute(n0));

      if(h0 < minHeight) minHeight = h0;


      Vector3 n1 = outer_product(points[2] - points[3], points[0] - points[3]);
      Scalar h1 = fabs(scalar_product(points[1] - points[3], n1) / absolute(n1));

      if(h1 < minHeight) minHeight = h1;

      Vector3 n2 = outer_product(points[3] - points[0], points[1] - points[0]);
      Scalar h2 = fabs(scalar_product(points[2] - points[0], n2) / absolute(n2));

      if(h2 < minHeight) minHeight = h2;

      Vector3 n3 = outer_product(points[0] - points[1], points[2] - points[1]);
      Scalar h3 = fabs(scalar_product(points[3] - points[1], n3) / absolute(n3));

      if(h3 < minHeight) minHeight = h3;
    }

    return minHeight;
  }

  Scalar ComputeAvancedTimeStep(ElasticProperties *submeshMaterialSettings, IndexType *submeshHeirarchialMultiplier)
  {
    Scalar minStep = Scalar(1e5);
    IndexType cells = this->GetCellsCount();
    for(IndexType i = 0; i < cells; i++)
    {
      IndexType submeshIndex = this->GetCellSubmeshIndex(i);

      typename GeomMesh::Cell *cell = this->GetCell(i);
      typename GeomMesh::Vector3 points[4];
      points[0] = this->GetNode(cell->incidentNodes[0])->pos;
      points[1] = this->GetNode(cell->incidentNodes[1])->pos;
      points[2] = this->GetNode(cell->incidentNodes[2])->pos;
      points[3] = this->GetNode(cell->incidentNodes[3])->pos;
      char tmpAxis;

      Scalar step = (ComputeInterpolationMinHeight(&tmpAxis, &(points[0].x)) / submeshMaterialSettings[submeshIndex].c1()) / submeshHeirarchialMultiplier[submeshIndex];

      if(step < minStep) minStep = step;
    }
    return minStep;
  }

  Scalar ComputeAvancedMinHeight()
  {
    Scalar minHeight = Scalar(1e5);
    IndexType cells = this->GetCellsCount();
    for(IndexType i = 0; i < cells; i++)
    {
      typename GeomMesh::Cell *cell = this->GetCell(i);
      typename GeomMesh::Vector3 points[4];
      points[0] = this->GetNode(cell->incidentNodes[0])->pos;
      points[1] = this->GetNode(cell->incidentNodes[1])->pos;
      points[2] = this->GetNode(cell->incidentNodes[2])->pos;
      points[3] = this->GetNode(cell->incidentNodes[3])->pos;
      char tmpAxis;

      Scalar height = ComputeInterpolationMinHeight(&tmpAxis, &(points[0].x));

      if(height < minHeight) minHeight = height;
    }
    return minHeight;
  }

  Scalar ComputeAdvancedMaxHeight()
  {
    Scalar maxHeight = 0;
    IndexType cells = this->GetCellsCount();
    for(IndexType i = 0; i < cells; i++)
    {
      typename GeomMesh::Cell *cell = this->GetCell(i);
      typename GeomMesh::Vector3 points[4];
      points[0] = this->GetNode(cell->incidentNodes[0])->pos;
      points[1] = this->GetNode(cell->incidentNodes[1])->pos;
      points[2] = this->GetNode(cell->incidentNodes[2])->pos;
      points[3] = this->GetNode(cell->incidentNodes[3])->pos;
      Scalar height = Interpolation_MaxAltitude_3D_Swift(cell->data.axis, InterpolationOrder, &(points[0].x));
      if(height > maxHeight) maxHeight = height;
    }
    return maxHeight;
  }
private:
  Scalar ComputeInterpolationMinHeight(char *bestAxis, Scalar *coords)
  {
    return Interpolation_MinAltitude_3D_Swift(bestAxis, InterpolationOrder, coords);
    //return Interpolation_MinAltitude_3D_Swift_interpolation_2_4(bestAxis, coords);
  }

  void ComputeCellAxes()
  {
    for(IndexType i = 0; i < this->GetCellsCount(); i++)
    {
      Vector3 points[4];
      for(IndexType j = 0; j < 4; j++)
      {
        points[j] = this->GetNode(this->GetCell(i)->incidentNodes[j])->pos;
      }

      ComputeInterpolationMinHeight(&(cellAxes[i]), &(points[0].x));
    };
  }
};
/*template<typename Scalar, typename Vector3, int InterpolationOrder, typename InterpolatorT>
template<>
void ElasticMeshT::LoadEdgesState<0>(IniStateMaker<Vector3, Elastic> *stateMaker){}*/

#include "ExternalInterpolator/interpolation.h"
#include <limits>


/*template<typename Scalar, typename Vector3>
bool PointInCell(const Vector3 points[4], Vector3 testPoint)
{
  Scalar eps = 0;//std::numeric_limits<float>::epsilon();//Scalar(1e-9);

  Scalar side0 = mixed_product(points[1] - points[0], points[2] - points[0], testPoint - points[0]);
  Scalar side1 = mixed_product(points[2] - points[0], points[3] - points[0], testPoint - points[0]);
  Scalar side2 = mixed_product(points[3] - points[0], points[1] - points[0], testPoint - points[0]);
  Scalar side3 = mixed_product(points[3] - points[1], points[2] - points[1], testPoint - points[1]);
  if (side0 >= -eps && side1 >= -eps && side2 >= -eps && side3 >= -eps) return 1;
  if (side0 <=  eps && side1 <=  eps && side2 <=  eps && side3 <=  eps) return 1;
  return 0;
}*/

/*template<typename Scalar, typename Vector3>
bool ProjectPointAgainstLine(const Vector3d &t1, const Vector3d &t2, const Vector3d &p, Vector3d &p0)
{
	Vector3d v1 = p - t1;
	Vector3d v2 = t2 - t1;
	p0 = t1 + v2 * (scalar_product(v1 * v2) / v2.SquareLen());
	if((v1 * v2 >= 0.0f) && ((v1 * v2) / (v2.SquareLen()) <= 1.0f))
	{
		return 1;
	}else
	{
		return 0;
	}
}

template<typename Scalar, typename Vector3>
bool PointProjectionInsideTriangle(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &point)
{
	Vector3 n = outer_product((p2 - p1), (p3 - p1));

	Vector3 n1 = outer_product((p2 - p1), n);
	Vector3 n2 = outer_product((p3 - p2), n);
	Vector3 n3 = outer_product((p1 - p3), n);

	Scalar proj1 = scalar_product((point - p2), n1);
	Scalar proj2 = scalar_product((point - p3), n2);
	Scalar proj3 = scalar_product((point - p1), n3);

	if(proj1 > 0)
		return 0;
	if(proj2 > 0)
		return 0;
	if(proj3 > 0)
		return 0;
	return 1;
}

template<typename Scalar, typename Vector3>
bool ProjectPointAgainstTriangle(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &point, Vector3 &projectionPoint)
{
	Vector3 n = outer_product((p2 - p1), (p3 - p1));
	Scalar mult = (scalar_product(p1, n) - scalar_product(point, n)) / scalar_product(n, n);
	projectionPoint = point + scalar_product(n, mult);

	return PointProjectionInsideTriangle(p1, p2, p3, projectionPoint);
};

static IndexType ProjectPointOnCell(Vector3 point, Vector3 *nodePositions, Vector3 &projectedPoint)
{
  if(!PointIsInCell(nodePositions[0], nodePositions[1], nodePositions[2], nodePositions[3], point)) return IndexType(-1);

  Scalar bestProjectionSqr = T(1e5);
  Vector3 bestProjectionPoint = Vector3(0, 0, 0);
  IndexType bestProjectionType = IndexType(-1);

  //*face projections*
  if(ProjectPointAgainstTriangle(nodePositions[0], nodePositions[1], nodePositions[2], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 0;
    }
  }
  if(ProjectPointAgainstTriangle(nodePositions[0], nodePositions[1], nodePositions[3], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 0;
    }
  }
  if(ProjectPointAgainstTriangle(nodePositions[1], nodePositions[2], nodePositions[3], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 0;
    }
  }
  if(ProjectPointAgainstTriangle(nodePositions[2], nodePositions[0], nodePositions[3], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 0;
    }
  }

  //edges
  if(ProjectPointAgainstLine(nodePositions[0], nodePositions[1], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 1;
    }
  }
  if(ProjectPointAgainstLine(nodePositions[1], nodePositions[2], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 1;
    }
  }
  if(ProjectPointAgainstLine(nodePositions[2], nodePositions[0], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 1;
    }
  }
  if(ProjectPointAgainstLine(nodePositions[0], nodePositions[3], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 1;
    }
  }
  if(ProjectPointAgainstLine(nodePositions[1], nodePositions[3], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 1;
    }
  }
  if(ProjectPointAgainstLine(nodePositions[2], nodePositions[3], point, testProjection)
  {
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = testProjection;
      bestProjectionType = 1;
    }
  }

  //nodes
  {
    testProjection = nodePositions[0];
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = node;
      bestProjectionType = 1;
    }
  }
  {
    testProjection = nodePositions[1];
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = node;
      bestProjectionType = 2;
    }
  }
  {
    testProjection = nodePositions[2];
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = node;
      bestProjectionType = 2;
    }
  }
  {
    testProjection = nodePositions[3];
    Scalar projectionSqr = sqr(testProjection - point);
    if(projectionSqr < bestProjectionSqr)
    {
      bestProjectionSqr = projectionSqr;
      bestProjectionPoint = node;
      bestProjectionType = 2;
    }
  }
  projectedPoint = bestProjectionPoint;
  return bestProjectionType;
}*/

template<typename GeomSpace, int InterpolationOrderT>
class InterpolatorT
{
public:
  typedef typename GeomSpace::Scalar                     Scalar;
  typedef typename GeomSpace::Vector3                    Vector3;
  typedef typename GeomSpace::IndexType                  IndexType;

  const static IndexType InterpolationOrder = InterpolationOrderT;

  static Scalar GetUnboundValue(Vector3 interpolationPoint, Vector3 *nodePositions, Scalar *interpolationData, char axis)
  {
    Vector3 points[5];
    points[0] = interpolationPoint;
    for(int i = 0; i < 4; i++)
    {
      points[i + 1] = nodePositions[i];
    }
    return Interpolate_3D<Scalar>(InterpolationOrderT, 'P', axis, &(points[0].x), interpolationData);
  }

  static Vector3 GetUnboundGradient(Vector3 interpolationPoint, Vector3 *nodePositions, Scalar *interpolationData, char axis)
  {
    Scalar step = T(1e-5);

    Vector3 gradient;
    gradient.x =  GetUnboundValue(interpolationPoint + Vector3( step, 0, 0), nodePositions, interpolationData, axis) -
                  GetUnboundValue(interpolationPoint + Vector3(-step, 0, 0), nodePositions, interpolationData, axis);
    gradient.y =  GetUnboundValue(interpolationPoint + Vector3(0,  step, 0), nodePositions, interpolationData, axis) -
                  GetUnboundValue(interpolationPoint + Vector3(0, -step, 0), nodePositions, interpolationData, axis);
    gradient.z =  GetUnboundValue(interpolationPoint + Vector3(0, 0,  step), nodePositions, interpolationData, axis) -
                  GetUnboundValue(interpolationPoint + Vector3(0, 0, -step), nodePositions, interpolationData, axis);
    gradient /= step;

    return gradient;
  }



  static Scalar GetPointQuality(Vector3 interpolationPoint, Vector3 *nodePositions, char axis)
  {
    Vector3 points[5];
    points[0] = interpolationPoint;
    for(int i = 0; i < 4; i++)
    {
      points[i + 1] = nodePositions[i];
    }

    Vector3 smallTetraPoints[4];
    Nodes_Coord_3D<Scalar>(&(smallTetraCoords[0].x), interplationOrderT, axis, nodePositions);

    assert(PointIsInCell(smallTetraCoords[0], smallTetraCoords[1], smallTetraCoords[2], smallTetraCoords[3], interpolationPoint));

    bool stillInside = 1;
    IndexType maxIterations = 10;

    Vector3 currPoint = interpolationPoint;
    for(IndexType iteration = 0; iteration < maxIterations; iteration++);
    {
      Vector3 gradient = GetUnboundGradient(currPoint, nodePositions, interpolationData, axis);
      Scalar value = GetUnboundValue(currPoint, nodePositions, interpolationData, axis);
      Scalar gradientLen = absolute(gradient);
      if(gradientLen < Scalar(1e-5))
      {
        if(fabs(value) < Scalar(1e-5))
          return Scalar(iteration);
        else
          return Scalar(1.0);
      }
      if(value > 0) 
        dir = Scalar(1.0);
      else
        dir = Scalar(-1.0);

      Scalar step = 1e-1f;
      currPoint += (gradient / gradientLen) * dir * step;

      if(ProjectPointOnCell(currPoint, smallTetraPoints, newPoint) == 3) //closest point is a node
      {
        return Scalar(iteration) / Scalar(maxIterations);
      }
    }
    return Scalar(1.0);
  }
  static Scalar Interpolate(Vector3 interpolationPoint, Vector3 *nodePositions, Scalar *interpolationData, char axis)
  {
    Vector3 points[5];
    points[0] = interpolationPoint;
    for(int i = 0; i < 4; i++)
    {
      points[i + 1] = nodePositions[i];
    }
//    return Interpolate_3D_Mono<Scalar>(axis, &(points[0].x), interpolationData);
//    return Interpolate_3D<Scalar>(InterpolationOrderT, 'L', axis, &(points[0].x), interpolationData);
    return Interpolate_3D<Scalar>(InterpolationOrderT, 'P', axis, &(points[0].x), interpolationData);
//    return Interpolate_3D_interpolation_2_4<Scalar>('M', axis, &(points[0].x), interpolationData);
   /* {
      Scalar unboundValue = Interpolate_3D<Scalar>(InterpolationOrderT, 'M', axis, &(points[0].x), interpolationData);
      Scalar linearValue  = Interpolate_3D<Scalar>(InterpolationOrderT, 'L', axis, &(points[0].x), interpolationData);
      Scalar coef = Scalar(0.0);
//      Scalar ratio = Scalar(1.0) / (Scalar(1.0) + coef * fabs(unboundValue - linearValue) / (fabs(unboundValue + linearValue) + 1e-5f));
      Scalar ratio = 0;//Scalar(1.0) - GetPointQuality(interpolationPoint, nodePositions, axis);
      //Scalar ratio = Scalar(1.0) / (Scalar(1.0) + coef * fabs(unboundValue - linearValue) / (fabs(unboundValue + linearValue) + 1e-5f));
      return unboundValue * ratio + linearValue * (Scalar(1.0) - ratio);
      //return unboundValue * ratio + linearValue * (Scalar(1.0) - ratio);
    }*/
    //return Interpolate_3D_interpolation_2_4<Scalar>('M', axis, &(points[0].x), interpolationData);
//    return Interpolate_3D_interpolation_2_4<Scalar>('P', '2', &(points[0].x), interpolationData);
  }

  static Vector3 GetEdgePointPos(Vector3 edgePoint1, Vector3 edgePoint2, IndexType pointNum)
  {
    Vector3 resultPoint;
    Vector3 edgePoints[2];

    edgePoints[0] = edgePoint1;
    edgePoints[1] = edgePoint2;

    Interpolation_3D_Get_Edge<Scalar>(&(resultPoint.x), (int)InterpolationOrder, &(edgePoints[0].x), (int)pointNum);

    return resultPoint;
  }

  static Vector3 GetFacePointPos(Vector3 facePoint1, Vector3 facePoint2, Vector3 facePoint3, IndexType pointNum)
  {
    Vector3 resultPoint;
    Vector3 facePoints[3];

    facePoints[0] = facePoint1;
    facePoints[1] = facePoint2;
    facePoints[2] = facePoint3;

    Interpolation_3D_Get_Side<Scalar>(&(resultPoint.x), (int)InterpolationOrder, &(facePoints[0].x), (int)pointNum);

    return resultPoint;
  }

  static Vector3 GetCellPointPos(Vector3 cellPoint1, Vector3 cellPoint2, Vector3 cellPoint3, Vector3 cellPoint4, IndexType pointNum)
  {
    Vector3 resultPoint;
    Vector3 cellPoints[4];

    cellPoints[0] = cellPoint1;
    cellPoints[1] = cellPoint2;
    cellPoints[2] = cellPoint3;
    cellPoints[3] = cellPoint4;

    Interpolation_3D_Get_Volume<Scalar>(&(resultPoint.x), (int)InterpolationOrder, &(cellPoints[0].x), (int)pointNum);

    return resultPoint;
  }

  static void RotateEdgeData(Scalar *edgeData, IndexType *targetIndices)
  {
    InterpolationPermitation_Edge(InterpolationOrderT, targetIndices, edgeData);
  }
};

template<typename FunctionSpace>
ElasticVolumeMesh<Space2, FunctionSpace>::ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  ElasticVolumeMeshCommon<Space2, FunctionSpace>(solver, tolerance, hierarchyLevelsCount)
{
}

template<typename FunctionSpace>
inline typename ElasticVolumeMesh<Space2, FunctionSpace>::Elastic 
ElasticVolumeMesh<Space2, FunctionSpace>::GetAverageEdgeElastic(IndexType cellIndex, IndexType edgeNumber) const
{
  return volumeMesh.GetEdgeAverageSolution(cellIndex, edgeNumber);
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::MakeSnapshot(Elastic* destData, 
  IndexType width, IndexType height, 
  Vector boxPoint1, Vector boxPoint2, bool halfStepSolution)
{
  for(IndexType i = 0; i < width * height; i++)
  {
    for(IndexType j = 0; j < volumeMesh.dimsCount; j++)
    {
      destData[i].values[j] = Scalar(0.0);
    }
  }

  for(IndexType i = 0; i < volumeMesh.cells.size(); i++)
  {
    Vector points[3];
    volumeMesh.GetCellVertices(i, points);

    Scalar 
      xmin = points[0].x, xmax = points[0].x,
      ymin = points[0].y, ymax = points[0].y;
    for(int j = 0; j < 3; j++)
    {
      if(points[j].x < xmin) xmin = points[j].x;
      if(points[j].x > xmax) xmax = points[j].x;
      if(points[j].y < ymin) ymin = points[j].y;
      if(points[j].y > ymax) ymax = points[j].y;
    }

    int ixmin = int(((xmin - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * Scalar(width  - 1));
    int ixmax = int(((xmax - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * Scalar(width  - 1)) + 1;
    int iymin = int(((ymin - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * Scalar(height - 1));
    int iymax = int(((ymax - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * Scalar(height - 1)) + 1;

    if(ixmin < 0) ixmin = 0;
    if(ixmin > int(width - 1)) ixmin = int(width - 1);
    if(ixmax < 0) ixmax = 0;
    if(ixmax > int(width - 1)) ixmax = int(width - 1);

    if(iymin < 0) iymin = 0;
    if(iymin > int(height - 1)) iymin = int(height - 1);
    if(iymax < 0) iymax = 0;
    if(iymax > int(height - 1)) iymax = int(height - 1);

    for(int x = ixmin; x <= ixmax; x++)
    {
      for(int y = iymin; y <= iymax; y++)
      {
        Vector point(boxPoint1.x + (boxPoint2.x - boxPoint1.x) * (Scalar(x) / Scalar(width - 1)),
                      boxPoint1.y + (boxPoint2.y - boxPoint1.y) * (Scalar(y) / Scalar(height - 1)));

        if (PointInCell<Scalar>(points, point))
        {
          Elastic e = InterpolateElastic(i, point, halfStepSolution);
          e.MakeDimension(tensionDimensionlessMult, velocityDimensionlessMult);
          destData[y * width + x] = e;
        }
      }
    }
  }
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::MakeSnapshot(Elastic* destData, 
  const Vector& origin, const Vector& spacing,
  const Vector& boxPoint1, const Vector& boxPoint2,
  bool halfStepSolution)
{
  IntAABB globalSnapshotArea;
  AABB snapshotArea;
  snapshotArea.Set(boxPoint1, boxPoint2);

  IndexVector snapshotSize = ::ComputeSnapshotAABBSize<Space2>(origin, spacing, snapshotArea, 
    &globalSnapshotArea);

  int globalIxmin = globalSnapshotArea.boxPoint1.x;
  int globalIymin = globalSnapshotArea.boxPoint1.y;

  for(IndexType elasticIndex = 0; elasticIndex < snapshotSize.GetVolume(); ++elasticIndex)
  {
    std::fill_n(destData[elasticIndex].values, volumeMesh.dimsCount, Scalar(0.0));
  }

  for(IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    Vector points[Space::NodesPerCell];
    volumeMesh.GetCellVertices(cellIndex, points);

    Scalar 
      xmin = points[0].x, xmax = points[0].x,
      ymin = points[0].y, ymax = points[0].y;

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if(points[nodeNumber].x < xmin) xmin = points[nodeNumber].x;
      if(points[nodeNumber].x > xmax) xmax = points[nodeNumber].x;
      if(points[nodeNumber].y < ymin) ymin = points[nodeNumber].y;
      if(points[nodeNumber].y > ymax) ymax = points[nodeNumber].y;
    }

    int ixmin = int((xmin - origin.x) / spacing.x) - 1;
    int ixmax = int((xmax - origin.x) / spacing.x) + 1;

    int iymin = int((ymin - origin.y) / spacing.y) - 1;
    int iymax = int((ymax - origin.y) / spacing.y) + 1;

    for(int x = ixmin; x <= ixmax; x++)
    {
      for(int y = iymin; y <= iymax; y++)
      {
        if(y - globalIymin >= 0 && y - globalIymin < int(snapshotSize.y) &&
           x - globalIxmin >= 0 && x - globalIxmin < int(snapshotSize.x))
        {
          Vector point(origin.x + spacing.x * Scalar(x),
                       origin.y + spacing.y * Scalar(y));

          if(PointInCell<Scalar>(points, point) && snapshotArea.Includes(point))
          {
            Elastic e = InterpolateElastic(cellIndex, point, halfStepSolution);
            e.MakeDimension(tensionDimensionlessMult, velocityDimensionlessMult);
            destData[(y - globalIymin) * snapshotSize.x + x - globalIxmin] = e;
          }
        }
      }
    }
  }
}

/*
  if a cell is destructed then only compress and shear tensions exist

  s:s < 2k^2 -- von Mises criterion;
  k = k0 - a * pressure, a > 0.
*/
template<typename MeshType>
struct ContinuousDestructionCorrector
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::Tensor    Tensor;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  ContinuousDestructionCorrector(MeshType* mesh, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->cellIndex = cellIndex;
  }

  void operator()(IndexType basisPointIndex, Scalar* values) const
  {
    bool brittle = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, basisPointIndex, MeshType::VolumeMeshTypeCommon::Basis);

    if (brittle)
    {
      const Scalar k0 = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.k0;
      const Scalar  a = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.a;
      const Scalar pressure = elastic.GetPressure();
      const Scalar k = k0 + a * pressure;

      if (mesh->ProcessPlasticity(k, elastic, false))
      {
        mesh->DestroyCell(cellIndex, mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
      }
    }

    if (mesh->volumeMesh.cellMediumParameters[cellIndex].destroyed)
    {
      Scalar mainStresses[2];
      Scalar maxTangentStress = Scalar(0.5) * (sqrt(Sqr(elastic.GetXX() - elastic.GetYY()) + 4 * Sqr(elastic.GetXY())));

      mainStresses[0] = Scalar(0.5) * (elastic.GetXX() + elastic.GetYY()) + maxTangentStress;
      mainStresses[1] = Scalar(0.5) * (elastic.GetXX() + elastic.GetYY()) - maxTangentStress;

      Vector normal = Vector(mainStresses[0] - elastic.GetYY(), elastic.GetXY()).GetNorm();

      for (IndexType stressIndex = 0; stressIndex < 2; ++stressIndex)
      {
        mainStresses[stressIndex] = std::min<Scalar>(0, mainStresses[stressIndex]);
      }

      // rotate from local main tensions axes to global coordinates
      elastic.SetTension(mainStresses[0] * Sqr(normal.x) + mainStresses[1] * Sqr(normal.y), 
        mainStresses[0] * Sqr(normal.y) + mainStresses[1] * Sqr(normal.x),
        (mainStresses[0] - mainStresses[1]) * normal.x * normal.y);
    }
    
    std::copy(elastic.values, elastic.values + mesh->volumeMesh.dimsCount, values);
  }

  MeshType* mesh;
  IndexType cellIndex;
};

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::HandleContinuousDestruction()
{
  typedef ContinuousDestructionCorrector< ElasticVolumeMesh<Space2, FunctionSpace> > DestructionCorrectorType;
  ApplyCorrector< ElasticVolumeMesh<Space2, FunctionSpace>, DestructionCorrectorType>(this);
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::UnfoldMesh(Scalar minHeight, IndexType iterationsCount)
{
  bool found = true;
  for(IndexType iterationIndex = 0; (iterationIndex < iterationsCount) && found; iterationIndex++)
  {
    found = false;
    for(IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); cellIndex++)
    {
      IndexType cellIndices[3];
      volumeMesh.GetFixedCellIndices(cellIndex, cellIndices);

      for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
      {
        Vector &currVertex = volumeMesh.nodes [cellIndices[nodeNumber]].pos;
        Vector &edgeVertex0 = volumeMesh.nodes[cellIndices[(nodeNumber + 1) % 3]].pos;
        Vector &edgeVertex1 = volumeMesh.nodes[cellIndices[(nodeNumber + 2) % 3]].pos;
        // Vector oppositeEdgeVertices[2];

        Scalar side;

        Vector projectionPoint;
        if(ProjectPointAgainstLine(edgeVertex0, edgeVertex1, currVertex, projectionPoint, side) &&
           (projectionPoint - currVertex).SquareLen() < Sqr(minHeight))
        {
          found = true;
          Scalar len0 = (projectionPoint - edgeVertex0).Len();
          Scalar len1 = (projectionPoint - edgeVertex1).Len();
          Scalar ratio = len0 / (len0 + len1);

          Scalar pushDepth = minHeight - (projectionPoint - currVertex).Len();
          Vector pushVector = (projectionPoint - currVertex).GetNorm() * pushDepth;

          Scalar coef = 0.5;
          currVertex  -= pushVector * 0.5 * coef;
          edgeVertex0 += pushVector * 0.5 * (Scalar(1.0) - ratio) * coef;
          edgeVertex1 += pushVector * 0.5 * (              ratio) * coef;
        }
      }
    }
  }
}

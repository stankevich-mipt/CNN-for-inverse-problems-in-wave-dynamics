template<typename FunctionSpace>
ElasticVolumeMesh<Space3, FunctionSpace>::ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  ElasticVolumeMeshCommon<Space3, FunctionSpace>(solver, tolerance, hierarchyLevelsCount)
{
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::UnfoldMesh(Scalar minHeight, IndexType iterationsCount)
{
}

template<typename MeshType>
struct ContinuousDestructionCorrector3
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::Tensor    Tensor;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  ContinuousDestructionCorrector3(MeshType* mesh, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->cellIndex = cellIndex;
  }

  void operator()(IndexType basisPointIndex, Scalar* values) const
  {
    bool brittle = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, basisPointIndex, MeshType::VolumeMeshType::Basis);

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
      Scalar principalStresses[Space3::Dimension];

      Scalar test[Space3::Dimension];
      elastic.GetTension().GetEigenValues(test);

      bool needToCorrect = false;
      for (IndexType stressIndex = 0; stressIndex < Space3::Dimension; ++stressIndex)
      {
        if (principalStresses[stressIndex] > 0)
        {
          needToCorrect = true;
          break;
        }
      }

      if (needToCorrect)
      {
        Vector principalNormals[Space3::Dimension];
        elastic.GetTension().GetEigenValues(principalStresses, principalNormals);

        for (IndexType stressIndex = 0; stressIndex < Space3::Dimension; ++stressIndex)
        {
          principalStresses[stressIndex] = std::min<Scalar>(0, principalStresses[stressIndex]);
        }

        Tensor t;
        t.xx = principalStresses[0];
        t.yy = principalStresses[1];
        t.zz = principalStresses[2];

        std::swap(principalNormals[0].y, principalNormals[1].x);
        std::swap(principalNormals[0].z, principalNormals[2].x);
        std::swap(principalNormals[1].z, principalNormals[2].y);

        t = t.RotateAxes(principalNormals[0], principalNormals[1], principalNormals[2]);
        // rotate from local principal tensions axes to global coordinates
        elastic.SetTension(t);
      }
    }

    std::copy(elastic.values, elastic.values + mesh->volumeMesh.dimsCount, values);
  }

  MeshType* mesh;
  IndexType cellIndex;
};

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::HandleContinuousDestruction()
{
  typedef ContinuousDestructionCorrector3< ElasticVolumeMesh<Space3, FunctionSpace> > DestructionCorrectorType;
  ApplyCorrector< ElasticVolumeMesh<Space3, FunctionSpace>, DestructionCorrectorType>(this);
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::MakeSnapshot(Elastic* destData,
  const Vector& origin, const Vector& spacing,
  const Vector& boxPoint1, const Vector& boxPoint2,
  bool halfStepSolution)
{
  IntAABB globalSnapshotArea;
  AABB snapshotArea;
  snapshotArea.Set(boxPoint1, boxPoint2);

  IndexVector snapshotSize = ::ComputeSnapshotAABBSize<Space3>(origin, spacing, snapshotArea,
    &globalSnapshotArea);

  int globalIxmin = globalSnapshotArea.boxPoint1.x;
  int globalIymin = globalSnapshotArea.boxPoint1.y;
  int globalIzmin = globalSnapshotArea.boxPoint1.z;

  for (IndexType elasticIndex = 0; elasticIndex < snapshotSize.GetVolume(); ++elasticIndex)
  {
    std::fill_n(destData[elasticIndex].values, volumeMesh.dimsCount, Scalar(0.0));
  }

  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    Vector points[Space::NodesPerCell];
    volumeMesh.GetCellVertices(cellIndex, points);

    Scalar
      xmin = points[0].x, xmax = points[0].x,
      ymin = points[0].y, ymax = points[0].y,
      zmin = points[0].z, zmax = points[0].z;

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (points[nodeNumber].x < xmin) xmin = points[nodeNumber].x;
      if (points[nodeNumber].x > xmax) xmax = points[nodeNumber].x;
      if (points[nodeNumber].y < ymin) ymin = points[nodeNumber].y;
      if (points[nodeNumber].y > ymax) ymax = points[nodeNumber].y;
      if (points[nodeNumber].z < zmin) zmin = points[nodeNumber].z;
      if (points[nodeNumber].z > zmax) zmax = points[nodeNumber].z;
    }

    int ixmin = int((xmin - origin.x) / spacing.x) - 1;
    int ixmax = int((xmax - origin.x) / spacing.x) + 1;

    int iymin = int((ymin - origin.y) / spacing.y) - 1;
    int iymax = int((ymax - origin.y) / spacing.y) + 1;

    int izmin = int((zmin - origin.z) / spacing.z) - 1;
    int izmax = int((zmax - origin.z) / spacing.z) + 1;

    for (int x = ixmin; x <= ixmax; x++)
    {
      for (int y = iymin; y <= iymax; y++)
      {
        for (int z = izmin; z <= izmax; z++)
        {
          if (x - globalIxmin >= 0 && x - globalIxmin < int(snapshotSize.x) &&
              y - globalIymin >= 0 && y - globalIymin < int(snapshotSize.y) &&              
              z - globalIzmin >= 0 && z - globalIzmin < int(snapshotSize.z))
          {
            Vector point(origin.x + spacing.x * Scalar(x),
                         origin.y + spacing.y * Scalar(y),
                         origin.z + spacing.z * Scalar(z));

            if (PointInCell<Scalar>(points, point) && snapshotArea.Includes(point))
            {
              Elastic e = InterpolateElastic(cellIndex, point, halfStepSolution);
              e.MakeDimension(tensionDimensionlessMult, velocityDimensionlessMult);
              destData[(z - globalIzmin) * snapshotSize.x * snapshotSize.y +
                       (y - globalIymin) * snapshotSize.x + 
                       (x - globalIxmin)] = e;
            }
          }
        }
      }
    }
  }
}

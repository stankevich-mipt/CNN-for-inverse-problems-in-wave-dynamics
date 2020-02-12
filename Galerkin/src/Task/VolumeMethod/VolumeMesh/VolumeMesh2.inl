template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::
  LoadGeom(Vector* vertexPositions, IndexType* cellIndices, IndexType verticesCount, IndexType cellsCount,
           EdgePairIndices* contactEdges, IndexType* contactEdgesCount, IndexType contactTypesCount,
           BoundaryEdge*    boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount,
           MediumParameters* mediumParameters,
           IndexType *internalContactTypes)
{
  printf("Building geom mesh topology \n");

  GeomMesh<Space>::LoadGeom(vertexPositions, cellIndices, verticesCount, cellsCount);
  GeomMesh<Space>::BuildTopologyInfos();
  GeomMesh<Space>::BuildAdditionalTopology(
    contactEdges,  contactEdgesCount,  contactTypesCount,
    boundaryEdges, boundaryEdgesCount, boundaryTypesCount,
    internalContactTypes);

  Initialize();
  printf("Building volume method additional matrices \n");
  BuildMatrices();

  // medium parameters setting
  if (mediumParameters) 
  {
    cellMediumParameters.resize(cellsCount);
    std::copy(mediumParameters, mediumParameters + cellsCount, cellMediumParameters.begin());
  }

  printf("Loading done \n");
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::BuildMatrices()
{
  #pragma omp parallel
  for(int functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
  {
    #pragma omp for
    for(int functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
    {
      printf(".");
      cellVolumeIntegrals(functionIndex0, functionIndex1) =
        functionSpace->ComputeCellVolumeIntegral       (functionIndex1, functionIndex0);

      Vector derivativeIntegral = functionSpace->ComputeDerivativeVolumeIntegral (functionIndex0, functionIndex1);
      xDerivativeVolumeIntegrals (functionIndex0, functionIndex1) = derivativeIntegral.x;
      yDerivativeVolumeIntegrals (functionIndex0, functionIndex1) = derivativeIntegral.y;

      for(IndexType srcEdgeNumber = 0; srcEdgeNumber < 3; srcEdgeNumber++)
      {
        outgoingFlux.srcEdges[srcEdgeNumber].surfaceIntegral(functionIndex0, functionIndex1) =
          functionSpace->ComputeOutgoingFlux(srcEdgeNumber, functionIndex0, functionIndex1);
        for(IndexType dstEdgeNumber = 0; dstEdgeNumber < 3; dstEdgeNumber++)
        {
          incomingFlux.srcEdges[srcEdgeNumber].dstEdges[dstEdgeNumber].surfaceIntegral(functionIndex0, functionIndex1) =
            functionSpace->ComputeIncomingFlux(srcEdgeNumber, dstEdgeNumber, functionIndex0, functionIndex1);
        }
      }
    }
  }

  cellVolumeIntegralsInv = cellVolumeIntegrals.inverse();

  xDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;
  yDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;

  #ifdef USE_SPARSE_MATRIX_FOR_DERIVATIVES
    xDerivativeVolumeIntegralsSparse.reserve(Eigen::VectorXi::Constant(functionsCount, functionsCount));
    yDerivativeVolumeIntegralsSparse.reserve(Eigen::VectorXi::Constant(functionsCount, functionsCount));

    for (int functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
    {
      for (int functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
      {
        if (fabs(xDerivativeVolumeIntegrals(functionIndex0, functionIndex1)) > std::numeric_limits<Scalar>::epsilon())
          xDerivativeVolumeIntegralsSparse.insert(functionIndex0, functionIndex1) = xDerivativeVolumeIntegrals(functionIndex0, functionIndex1);
        if (fabs(yDerivativeVolumeIntegrals(functionIndex0, functionIndex1)) > std::numeric_limits<Scalar>::epsilon())
          yDerivativeVolumeIntegralsSparse.insert(functionIndex0, functionIndex1) = yDerivativeVolumeIntegrals(functionIndex0, functionIndex1);
      }
    }

    xDerivativeVolumeIntegralsSparse.makeCompressed();
    yDerivativeVolumeIntegralsSparse.makeCompressed();
  #endif

  for(IndexType srcEdgeNumber = 0; srcEdgeNumber < 3; srcEdgeNumber++)
  {
    outgoingFlux.srcEdges[srcEdgeNumber].surfaceIntegral *= cellVolumeIntegralsInv;

    for(IndexType dstEdgeNumber = 0; dstEdgeNumber < 3; dstEdgeNumber++)
    {
      incomingFlux.srcEdges[srcEdgeNumber].dstEdges[dstEdgeNumber].surfaceIntegral *= cellVolumeIntegralsInv;
    }
  }

  Eigen::Matrix<Scalar, functionsCount, functionsCount> testMatrix;
  Scalar err = 0;

  testMatrix = cellVolumeIntegrals * cellVolumeIntegralsInv;
  for(IndexType i = 0; i < functionsCount; i++)
  {
    for(IndexType j = 0; j < functionsCount; j++)
    {
      if(i == j)
        err += fabs(testMatrix(i, j) - Scalar(1.0));
      else
        err += fabs(testMatrix(i, j));
    }
  }
  printf("\nPrecomputations complete, volume integral matrix error is : %f\n", err);

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    cellVolumeAverageIntegrals[functionIndex] = functionSpace->ComputeCellVolumeIntegral(functionIndex);

    for(IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
    {
      edgeAverages[edgeNumber].surfaceIntegral[functionIndex] =
        functionSpace->ComputeEdgeFlux(edgeNumber, functionIndex);
    }
  }
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::
  GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState)
{
  #pragma omp parallel 
  {
    int threadIndex  = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    MatrixXDimFunc currCellValues;
    MatrixXDimFunc correspondingCellValues;
    MatrixXDimFunc timeDerivatives;
    MatrixXDimFunc edgeFlux;

    MatrixXDim edgeTransformMatrix;
    MatrixXDim edgeTransformMatrixInv;

    MatrixXDim xnAuxMatrix;
    MatrixXDim xInteriorMatrix;
    MatrixXDim xExteriorMatrix;

    MatrixXDim xMatrix;
    MatrixXDim yMatrix;

    MatrixXDim xMixedMatrix;
    MatrixXDim yMixedMatrix;

    Eigen::Matrix<Scalar, 1, dimsCount> boundaryMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> leftContactMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> rightContactMatrix;

    MatrixXDimFunc boundaryInfoValues;
    MatrixXDimFunc sourceValues;
    MatrixXDimFunc sourcePointValues;

    MatrixXDimFunc flux;

    Vector cellVertices[Space::NodesPerCell];
    Scalar tmp[dimsCount];

    for(int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      /*
        There are regular cells, where we compute solution,
        and domain boundary cells which are taken from neighbouring domains and should not be computed here.
      */

      bool auxCell;
      if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell)) continue;

      timeDerivatives.setZero();

      if (IsCellRegular(cellIndex) && isCellAvailable[cellIndex] && !cellMediumParameters[cellIndex].fixed)
      {
        bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
        for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
        {
          for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
          {
            currCellValues(valueIndex, functionIndex) =
              useHalfStepSolution ?
              halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
              cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];
          }
        }

        GetCellVertices(cellIndex, cellVertices);
        Scalar invJacobian = Scalar(1.0) / fabs(GetCellDeformJacobian(cellVertices));

        IndexType cellIncidentNodes[Space::NodesPerCell];
        GetFixedCellIndices(cellIndex, cellIncidentNodes);

        system.BuildXMatrix(cellMediumParameters[cellIndex], xMatrix);
        system.BuildYMatrix(cellMediumParameters[cellIndex], yMatrix);

        for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
        {
          edgeFlux.setZero();

          IndexType edgeNodeIndices[Space::NodesPerEdge];
          GetCellEdgeNodes(cellIncidentNodes, edgeNumber, edgeNodeIndices);

          Vector edgeGlobalVertices[Space::NodesPerEdge];
          for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; nodeNumber++)
          {
            edgeGlobalVertices[nodeNumber] = nodes[edgeNodeIndices[nodeNumber]].pos;
          }

          system.BuildEdgeTransformMatrix(edgeGlobalVertices, edgeTransformMatrix);
          system.BuildEdgeTransformMatrixInv(edgeGlobalVertices, edgeTransformMatrixInv);

          Scalar edgeLen = (edgeGlobalVertices[1] - edgeGlobalVertices[0]).Len(); //GetFaceSquare(faceGlobalVertices);
          Vector edgeNormal = GetEdgeExternalNormal(cellIndex, edgeNumber);

          IndexType correspondingCellIndex = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex;
          IndexType correspondingEdgeNumber = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber;
          IndexType interactionType = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;

          if (interactionType == IndexType(-1)) continue;

          system.BuildXnAuxMatrix(cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            xnAuxMatrix);

          // interior matrix
          system.BuildXnInteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            edgeNormal, xnAuxMatrix, xInteriorMatrix);

          // exterior matrix
          system.BuildXnExteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex], 
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            edgeNormal, xnAuxMatrix, xExteriorMatrix);

          if (correspondingCellIndex == IndexType(-1))
          {
            // bounary condition
            IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);

            // external force/velocity 
            BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);
            if (functor)
            {
              typedef BoundaryFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > FunctorWrapper;
              FunctorWrapper wrapper(functor, time, this, cellIndex, edgeNumber);

              functionSpace->template Decompose< FunctorWrapper, dimsCount >(wrapper, boundaryInfoValues.data());

              edgeFlux.noalias() += Scalar(2.0) * (xExteriorMatrix * edgeTransformMatrixInv * boundaryInfoValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral);
            }

            if (!allowDynamicCollisions || dynamicContactType == IndexType(-1)) //set 1 for regular boundary, 0 for dynamic collisions
            {
              system.BuildBoundaryMatrix(interactionType, boundaryMatrix);
              edgeFlux.noalias() +=  (xInteriorMatrix + xExteriorMatrix * boundaryMatrix.asDiagonal()) * 
                edgeTransformMatrixInv * currCellValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral;
            }
            else
            {
              GhostCellFunctionGetter<VolumeMeshT> functionGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::Solution);

              GhostCellFunctionGetter<VolumeMeshT> paramsGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::MediumParams);

              flux.setZero();

              // quadrature integration of numerical flux
              for (IndexType pointIndex = 0; pointIndex < quadraturePointsForBorder.size(); ++pointIndex)
              {
                Vector globalPoint = (edgeGlobalVertices[1] - edgeGlobalVertices[0]) * quadraturePointsForBorder[pointIndex] + edgeGlobalVertices[0];

                Vector refPoint = GlobalToRefVolumeCoords(globalPoint, cellVertices);
                typename System::ValueType interiorSolution = GetRefCellSolution(cellIndex, refPoint);
                
                MatrixMulVector(edgeTransformMatrixInv.data(), interiorSolution.values, tmp, dimsCount, dimsCount);
                std::copy(tmp, tmp + dimsCount, interiorSolution.values);

                MediumParameters exteriorParams;
                std::fill(exteriorParams.params, exteriorParams.params + MediumParameters::ParamsCount, 0);
                IndexType collidedCellIndex = IndexType(-1);

                if (collisionWidth > std::numeric_limits<Scalar>::epsilon())
                  collidedCellIndex = paramsGetter(globalPoint, exteriorParams.params);
                else
                {
                  if (collisionsInfo.collisionNodes[cellIndex].count > 0)
                  {
                    collidedCellIndex = paramsGetter(globalPoint,
                      collisionsInfo.pool.data() + collisionsInfo.collisionNodes[cellIndex].offset,
                      collisionsInfo.collisionNodes[cellIndex].count, exteriorParams.params);
                  }
                }

                typename System::ValueType exteriorSolution;

                if (collidedCellIndex != IndexType(-1))
                {
                  if (!functionGetter.TryGhostCell(globalPoint + edgeNormal * collisionWidth, collidedCellIndex, exteriorSolution.values))
                  {
                    assert(0);
                  }

                  MatrixMulVector(edgeTransformMatrixInv.data(), exteriorSolution.values, tmp, dimsCount, dimsCount);
                  std::copy(tmp, tmp + dimsCount, exteriorSolution.values);
                }

                typename System::ValueType riemannSolution =
                  system.GetRiemannSolution(interiorSolution, exteriorSolution,
                    cellMediumParameters[cellIndex], exteriorParams,
                    interactionType, // boundary type
                    dynamicContactType);

                for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
                {
                  Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refPoint, functionIndex);
                  for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
                  {
                    flux(valueIndex, functionIndex) +=
                      quadratureWeightsForBorder[pointIndex] *
                      riemannSolution.values[valueIndex] * basisFunctionValue;
                  }
                }
              }

              edgeFlux.noalias() += xMatrix * flux * cellVolumeIntegralsInv;
            }
          }
          else
          {
            // contact condition
            bool useHalfStepSolutionForCorrespondingCell = timeHierarchyLevelsManager.UseHalfStepSolutionForNeighbour(
              cellIndex, solverState, auxCell, correspondingCellIndex);

            for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
            {
              for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
              {
                correspondingCellValues(valueIndex, functionIndex) =
                  useHalfStepSolutionForCorrespondingCell ?
                  halfStepCellSolutions[correspondingCellIndex].basisVectors[functionIndex].values[valueIndex] :
                  cellSolutions[correspondingCellIndex].basisVectors[functionIndex].values[valueIndex];
              }
            }

            system.BuildContactMatrices(interactionType, leftContactMatrix, rightContactMatrix);

            // for glue contact left matrix equals 0
            if (!leftContactMatrix.isZero(std::numeric_limits<Scalar>::epsilon()))
            {
              edgeFlux.noalias() += xExteriorMatrix * leftContactMatrix.asDiagonal() * edgeTransformMatrixInv *
                currCellValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral;
            }

            // exterior side contribution
            edgeFlux.noalias() += xExteriorMatrix * rightContactMatrix.asDiagonal() * edgeTransformMatrixInv *
              correspondingCellValues * incomingFlux.srcEdges[edgeNumber].dstEdges[correspondingEdgeNumber].surfaceIntegral;

            // outgoing flux
            edgeFlux.noalias() += xInteriorMatrix * edgeTransformMatrixInv *
              currCellValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral;
          }
          timeDerivatives.noalias() += (edgeLen * invJacobian) * (edgeTransformMatrix * edgeFlux);
        }

        Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
        Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;

        xMixedMatrix.noalias() = xMatrix * refXDerivatives.x + yMatrix * refXDerivatives.y;
        yMixedMatrix.noalias() = xMatrix * refYDerivatives.x + yMatrix * refYDerivatives.y;

        #ifdef USE_SPARSE_MATRIX_FOR_DERIVATIVES
        timeDerivatives.noalias() -=
          xMixedMatrix * currCellValues * xDerivativeVolumeIntegralsSparse +
          yMixedMatrix * currCellValues * yDerivativeVolumeIntegralsSparse;
        #else
        timeDerivatives.noalias() -=
          xMixedMatrix * currCellValues * xDerivativeVolumeIntegrals +
          yMixedMatrix * currCellValues * yDerivativeVolumeIntegrals;
        #endif

        typename SystemT::SourceFunctorT* sourceFunctor = system.GetSourceFunctor();
        if (sourceFunctor)
        {
          typedef SourceFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > SourceFunctorWrapper;
          SourceFunctorWrapper wrapper(sourceFunctor, time, this, cellVertices);
          functionSpace->template Decompose< SourceFunctorWrapper, dimsCount >(wrapper, sourceValues.data());
          timeDerivatives.noalias() -= sourceValues;
        }

        // point sources
        for (IndexType sourceIndex = 0; sourceIndex < system.pointSources.size(); ++sourceIndex)
        {
          Vector point = system.pointSources[sourceIndex]->GetPoint();
          if (PointInCell<Scalar>(cellVertices, point))
          {
            Vector refPoint = GlobalToRefVolumeCoords(point, cellVertices);
            Scalar values[dimsCount];
            (*system.pointSources[sourceIndex])(time, values);
            for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
            {
              for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
              {
                sourcePointValues(valueIndex, functionIndex) =
                  values[valueIndex] * functionSpace->GetBasisFunctionValue(refPoint, functionIndex);
              }
            }
            timeDerivatives.noalias() -= sourcePointValues * cellVolumeIntegralsInv;
          }
        }
      }

      for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          derivatives[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex]
            = -timeDerivatives(valueIndex, functionIndex);
        }
      }
      ++targetCellIndex;
    }
  }
}

template<typename FunctionSpace, typename System>
typename Space2::Scalar VolumeMesh<Space2, FunctionSpace, System>::
  GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const
{
  return (cellVertices[1].x - cellVertices[0].x) * (cellVertices[2].y - cellVertices[0].y) -
         (cellVertices[2].x - cellVertices[0].x) * (cellVertices[1].y - cellVertices[0].y);
}

template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GetRefXDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dξ/dx, dξ/dy, dξ/dz) * J
{
    return Vector(
      cellVertices[2].y - cellVertices[0].y,
      cellVertices[0].x - cellVertices[2].x);
}

template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GetRefYDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dη/dx, dη/dy, dη/dz) * J
{
    return Vector(
      cellVertices[0].y - cellVertices[1].y,
      cellVertices[1].x - cellVertices[0].x);
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::GetRefDerivatives(
  Vector cellVertices[Space::NodesPerCell],
  Vector* refDerivatives) const
{
  Scalar invJacobian = Scalar(1.0) / fabs(GetCellDeformJacobian(cellVertices));
  refDerivatives[0] = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
  refDerivatives[1] = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;
}


template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const //x -> ?
{
  Scalar invJacobian = Scalar(1.0) / GetCellDeformJacobian(cellVertices);
  Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices);
  Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices);

  return invJacobian * Vector(
    cellVertices[2].x * cellVertices[0].y - cellVertices[0].x * cellVertices[2].y +
    globalCoords.x * refXDerivatives.x + globalCoords.y * refXDerivatives.y,

    cellVertices[0].x * cellVertices[1].y - cellVertices[1].x * cellVertices[0].y +
    globalCoords.x * refYDerivatives.x + globalCoords.y * refYDerivatives.y);
}

template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const // ? -> x
{
  return Vector(
    cellVertices[0].x + (cellVertices[1].x - cellVertices[0].x) * refCoords.x +
                        (cellVertices[2].x - cellVertices[0].x) * refCoords.y,
    cellVertices[0].y + (cellVertices[1].y - cellVertices[0].y) * refCoords.x +
                        (cellVertices[2].y - cellVertices[0].y) * refCoords.y
    );
}

template<typename FunctionSpace, typename System>
typename System::ValueType VolumeMesh<Space2, FunctionSpace, System>::
  GetEdgeAverageSolution(IndexType cellIndex, IndexType edgeNumber) const
{
  typename System::ValueType result(Scalar(0.0));

  for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
  {
    for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++) 
    {
      Scalar basisFunctionCoefficient = 
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];

      result.values[valueIndex] += basisFunctionCoefficient * edgeAverages[edgeNumber].surfaceIntegral[functionIndex];
    }
  }

  return result;
}

template<typename FunctionSpace, typename System>
typename System::ValueType VolumeMesh<Space2, FunctionSpace, System>::
GetFaceAverageSolution(IndexType cellIndex, IndexType faceNumber) const
{
  return GetEdgeAverageSolution(cellIndex, faceNumber);
}

template<typename FunctionSpace, typename System>
bool VolumeMesh<Space2, FunctionSpace, System>::IsCellRegular(IndexType cellIndex) const
{
  bool regularCell = true;
  for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
  {
    IndexType interactionType = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;
    if (interactionType == IndexType(-1)) regularCell = false;
  }
  return regularCell;
}

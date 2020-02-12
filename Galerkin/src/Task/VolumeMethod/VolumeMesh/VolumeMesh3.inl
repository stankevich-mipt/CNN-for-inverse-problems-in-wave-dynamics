#include "../../../Maths/QuadraturePrecomputer.h"

template<typename FunctionSpace, typename System>
void VolumeMesh<Space3, FunctionSpace, System>::
  LoadGeom(Vector* vertexPositions, IndexType* cellIndices, IndexType verticesCount, IndexType cellsCount,
           FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
           BoundaryFace*    boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount,
           MediumParameters* mediumParameters,
           IndexType *internalContactTypes)
{
  printf("Building geom mesh topology \n");
  GeomMesh<Space>::LoadGeom(vertexPositions, cellIndices, verticesCount, cellsCount);
  GeomMesh<Space>::BuildTopologyInfos();
  GeomMesh<Space>::BuildAdditionalTopology(
    contactFaces, contactFacesCount, contactTypesCount,
    boundaryFaces, boundaryFacesCount, boundaryTypesCount, 
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
void VolumeMesh<Space3, FunctionSpace, System>::BuildMatrices()
{
  #ifdef USE_DYNAMIC_MATRICIES
  xDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  yDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);

  for (IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
  {
    outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral.resize(functionsCount, functionsCount);
    for (IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
    {
      for (IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
      {
        incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral.resize(functionsCount, functionsCount);
      }
    }
  }

  xDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  yDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  zDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  #endif

  #pragma omp parallel for
  for(int functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
  {
    for(int functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
    {
      printf(".");
      fflush (stdout);
      cellVolumeIntegrals(functionIndex0, functionIndex1) =
        functionSpace->ComputeCellVolumeIntegral(functionIndex1, functionIndex0);

      Vector derivativeIntegral = functionSpace->ComputeDerivativeVolumeIntegral(functionIndex0, functionIndex1);
      xDerivativeVolumeIntegrals(functionIndex0, functionIndex1) = derivativeIntegral.x;
      yDerivativeVolumeIntegrals(functionIndex0, functionIndex1) = derivativeIntegral.y;
      zDerivativeVolumeIntegrals(functionIndex0, functionIndex1) = derivativeIntegral.z;

      for(IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
      {
        outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral(functionIndex0, functionIndex1) =
          functionSpace->ComputeOutgoingFlux(srcFaceNumber, functionIndex0, functionIndex1);
        for(IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
        {
          for(IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
          {
            incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral(functionIndex0, functionIndex1) =
              functionSpace->ComputeIncomingFlux(srcFaceNumber, dstFaceNumber, orientationNumber, functionIndex0, functionIndex1);
          }
        }
      }
    }
  }
  cellVolumeIntegralsInv = cellVolumeIntegrals.inverse();

  xDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;
  yDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;
  zDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;

  for(IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
  {
    outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral *= cellVolumeIntegralsInv;
    for(IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
    {
      for(IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
      {
        incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral *= cellVolumeIntegralsInv;
      }
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
    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
    {
      faceAverages[faceNumber].surfaceIntegral[functionIndex] =
        functionSpace->ComputeFaceFlux(faceNumber, functionIndex);
    }
  }
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space3, FunctionSpace, System>::
GetCurrDerivatives(Scalar *derivatives, const SolverState& solverState)
{
  Scalar res[32] = { 0 };

  #pragma omp parallel 
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    MatrixXDimFunc currCellValues;
    MatrixXDimFunc correspondingCellValues;
    MatrixXDimFunc timeDerivatives;
    MatrixXDimFunc faceFlux;

    MatrixXDim faceTransformMatrix;
    MatrixXDim faceTransformMatrixInv;

    MatrixXDim xnAuxMatrix;
    MatrixXDim xInteriorMatrix;
    MatrixXDim xExteriorMatrix;

    Eigen::Matrix<Scalar, 1, dimsCount> boundaryMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> leftContactMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> rightContactMatrix;

    MatrixXDim xMatrix;
    MatrixXDim yMatrix;
    MatrixXDim zMatrix;

    MatrixXDim xMixedMatrix;
    MatrixXDim yMixedMatrix;
    MatrixXDim zMixedMatrix;

    MatrixXDimFunc boundaryInfoValues;
    MatrixXDimFunc sourceValues;
    MatrixXDimFunc sourcePointValues;

    Scalar tmp[dimsCount];
    MatrixXDimFunc flux;

    #ifdef USE_DYNAMIC_MATRICIES
    currCellValues.resize(dimsCount, functionsCount);
    correspondingCellValues.resize(dimsCount, functionsCount);
    timeDerivatives.resize(dimsCount, functionsCount);

    faceTransformMatrix.resize(dimsCount, dimsCount);
    faceTransformMatrixInv.resize(dimsCount, dimsCount);

    xnAuxMatrix.resize(dimsCount, dimsCount);
    xInteriorMatrix.resize(dimsCount, dimsCount);
    xExteriorMatrix.resize(dimsCount, dimsCount);

    xMatrix.resize(dimsCount, dimsCount);
    yMatrix.resize(dimsCount, dimsCount);
    zMatrix.resize(dimsCount, dimsCount);

    xMixedMatrix.resize(dimsCount, dimsCount);
    yMixedMatrix.resize(dimsCount, dimsCount);
    zMixedMatrix.resize(dimsCount, dimsCount);

    boundaryInfoValues.resize(dimsCount, functionsCount);
    ghostValues.resize(dimsCount, functionsCount);
    sourceValues.resize(dimsCount, functionsCount);
    sourcePointValues.resize(dimsCount, functionsCount);

    flux.resize(dimsCount, functionsCount);
    #endif

    Vector cellVertices[Space::NodesPerCell];

    double beginPhase = MPI_Wtime();
    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
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

        system.BuildXMatrix(cellMediumParameters[cellIndex], xMatrix);
        system.BuildYMatrix(cellMediumParameters[cellIndex], yMatrix);
        system.BuildZMatrix(cellMediumParameters[cellIndex], zMatrix);

        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
        {
          Vector faceGlobalVertices[Space::NodesPerFace];
          GetCellFaceVertices(cellIndex, faceNumber, faceGlobalVertices);

          system.BuildFaceTransformMatrix(faceGlobalVertices, faceTransformMatrix);
          system.BuildFaceTransformMatrixInv(faceGlobalVertices, faceTransformMatrixInv);

          Scalar faceDeformJacobian = Scalar(2.0) * GetFaceSquare(faceGlobalVertices);

          Vector faceNormal = GetFaceExternalNormal(faceGlobalVertices);

          IndexType correspondingCellIndex = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
          IndexType correspondingFaceNumber = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber;
          IndexType correspondingFaceOrientation = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].orientation;
          IndexType interactionType = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;

          if (interactionType == IndexType(-1)) continue;

          faceFlux.setZero();

          system.BuildXnAuxMatrix(cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            xnAuxMatrix);

          // interior matrix
          system.BuildXnInteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xnAuxMatrix, xInteriorMatrix);

          // exterior matrix
          system.BuildXnExteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex], 
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xnAuxMatrix, xExteriorMatrix);

          if (correspondingCellIndex == IndexType(-1))
          {
            // bounary condition
            IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);

            // external force/velocity 
            BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);
            if (functor)
            {
              typedef BoundaryFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > FunctorWrapper;
              FunctorWrapper wrapper(functor, time, this, cellIndex, faceNumber);

              // TODO: replace it with integration over face
              functionSpace->template Decompose< FunctorWrapper, dimsCount >(wrapper, boundaryInfoValues.data());
              faceFlux.noalias() += Scalar(2.0) * xExteriorMatrix * faceTransformMatrixInv * boundaryInfoValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
            }

            if (!allowDynamicCollisions || dynamicContactType == IndexType(-1))
            {
              system.BuildBoundaryMatrix(interactionType, boundaryMatrix);
              // regular boundary          
              faceFlux.noalias() += (xInteriorMatrix + xExteriorMatrix * boundaryMatrix.asDiagonal()) *
                faceTransformMatrixInv * currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
            }
            else
            {
              // dynamic collision
              GhostCellFunctionGetter<VolumeMeshT> functionGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::Solution);

              // functionSpace->template Decompose< GhostCellFunctionGetter<VolumeMeshT>, dimsCount >(functionGetter, ghostValues.data());

              GhostCellFunctionGetter<VolumeMeshT> paramsGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::MediumParams);

              flux.setZero();

              // quadrature integration of numerical flux
              for (IndexType pointIndex = 0; pointIndex < quadraturePointsForBorder.size(); ++pointIndex)
              {
                Vector globalPoint = (faceGlobalVertices[1] - faceGlobalVertices[0]) * quadraturePointsForBorder[pointIndex].x +
                  (faceGlobalVertices[2] - faceGlobalVertices[0]) * quadraturePointsForBorder[pointIndex].y +
                  faceGlobalVertices[0];

                Vector refPoint = GlobalToRefVolumeCoords(globalPoint, cellVertices);
                typename System::ValueType interiorSolution = GetRefCellSolution(cellIndex, refPoint);

                MatrixMulVector(faceTransformMatrixInv.data(), interiorSolution.values, tmp, dimsCount, dimsCount);
                std::copy(tmp, tmp + dimsCount, interiorSolution.values);

                MediumParameters exteriorParams;
                std::fill(exteriorParams.params, exteriorParams.params + MediumParameters::ParamsCount, 0);
                IndexType collidedCellIndex = IndexType(-1);

                if (collisionWidth > std::numeric_limits<Scalar>::epsilon())
                  collidedCellIndex = paramsGetter(globalPoint + faceNormal * collisionWidth, exteriorParams.params);
                else
                {
                  if (collisionsInfo.collisionNodes[cellIndex].count > 0)
                  {
                    collidedCellIndex = paramsGetter(globalPoint + faceNormal * collisionWidth,
                      collisionsInfo.pool.data() + collisionsInfo.collisionNodes[cellIndex].offset,
                      collisionsInfo.collisionNodes[cellIndex].count, exteriorParams.params);
                  }
                }

                typename System::ValueType exteriorSolution; //GetRefCellSolution(ghostValues.data(), ghostRefPoint);
                if (collidedCellIndex != IndexType(-1))
                {
                  if (!functionGetter.TryGhostCell(globalPoint + faceNormal * collisionWidth, collidedCellIndex, exteriorSolution.values))
                  {
                    assert(0);
                  }

                  MatrixMulVector(faceTransformMatrixInv.data(), exteriorSolution.values, tmp, dimsCount, dimsCount);
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

              faceFlux.noalias() += xMatrix * flux * cellVolumeIntegralsInv;
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

            /*
            // README: uncomment when glide internal contact is needed
            system.BuildContactMatrices(interactionType, leftContactMatrix, rightContactMatrix);
            // for glue contact left matrix equals 0
            if (!leftContactMatrix.isZero(std::numeric_limits<Scalar>::epsilon()))
            {
            // interior side contribution
            faceFlux.noalias() += xExteriorMatrix * leftContactMatrix.asDiagonal() *
            faceTransformMatrixInv * currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
            }
            */

            // exterior side contribution
            faceFlux.noalias() += xExteriorMatrix * // rightContactMatrix.asDiagonal() *
              faceTransformMatrixInv * correspondingCellValues *
              incomingFlux.srcFaces[faceNumber].dstFaces[correspondingFaceNumber].orientations[correspondingFaceOrientation].surfaceIntegral;

            // outgoing flux
            faceFlux.noalias() += xInteriorMatrix * faceTransformMatrixInv *
              currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
          }

          timeDerivatives.noalias() += (faceDeformJacobian * invJacobian) * (faceTransformMatrix * faceFlux);
        }

        Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
        Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;
        Vector refZDerivatives = GetRefZDerivativesMulJacobian(cellVertices) * invJacobian;

        xMixedMatrix.noalias() = xMatrix * refXDerivatives.x + yMatrix * refXDerivatives.y + zMatrix * refXDerivatives.z;
        yMixedMatrix.noalias() = xMatrix * refYDerivatives.x + yMatrix * refYDerivatives.y + zMatrix * refYDerivatives.z;
        zMixedMatrix.noalias() = xMatrix * refZDerivatives.x + yMatrix * refZDerivatives.y + zMatrix * refZDerivatives.z;

        timeDerivatives.noalias() -=
          xMixedMatrix * currCellValues * xDerivativeVolumeIntegrals +
          yMixedMatrix * currCellValues * yDerivativeVolumeIntegrals +
          zMixedMatrix * currCellValues * zDerivativeVolumeIntegrals;

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

      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          derivatives[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex]
            = -timeDerivatives(valueIndex, functionIndex);
        }
      }
      ++targetCellIndex;
    }
    double endPhase = MPI_Wtime();
    res[threadIndex] = endPhase - beginPhase;
  }

  /*
  Scalar minTime = res[0];
  Scalar maxTime = res[0];
  int minIndex = 0;
  const int threadsCount = 8;
  for (int i = 0; i < threadsCount; ++i)
  {
    if (res[i] < minTime)
    {
      minTime = res[i];
      minIndex = i;
    }
    maxTime = std::max(res[i], maxTime);

  }
  std::cout << "\n";

  if (solverState.globalStepIndex % 20 == 0)
  {
    std::cout << "VolumeMesh: GetCurrDerivatives " << (maxTime - minTime) / minTime << " " << minIndex << "\n";

    for (IndexType threadIndex = 0; threadIndex < threadsCount; ++threadIndex)
    {
      IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
      int segmentBegin = threadSegmentBegins[stateIndex];
      int segmentEnd = threadSegmentEnds[stateIndex];

      int count = 0;

      for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
      {
        bool auxCell;
        if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell)) continue;

        if (IsCellRegular(cellIndex) && isCellAvailable[cellIndex] && !cellMediumParameters[cellIndex].fixed)
        {
          //count += 6 + 4;
          count++;
        }
        continue;

        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
        {
          IndexType correspondingCellIndex = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
          if (correspondingCellIndex != IndexType(-1))
          {
            count += 6;
          }
          else
          {
            count += 3;
          }
        }

      }
      std::cout << count << " ";
    }
    std::cout << "\n";

    for (IndexType threadIndex = 0; threadIndex < 8; ++threadIndex)
    {
      std::cout << res[threadIndex] << " ";
    }
    std::cout << "\n";
  }
  */
}

template< typename FunctionSpace, typename System>
typename Space3::Scalar VolumeMesh<Space3, FunctionSpace, System>::
GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const
{
  // deform jacobian equals 6 * tetrahedron volume
  return
    cellVertices[0].x * (cellVertices[1].y * (cellVertices[3].z - cellVertices[2].z) +
    cellVertices[2].y * (cellVertices[1].z - cellVertices[3].z) +
    cellVertices[3].y * (cellVertices[2].z - cellVertices[1].z)) +

    cellVertices[1].x * (cellVertices[0].y * (cellVertices[2].z - cellVertices[3].z) +
    cellVertices[2].y * (cellVertices[3].z - cellVertices[0].z) +
    cellVertices[3].y * (cellVertices[0].z - cellVertices[2].z)) +

    cellVertices[2].x * (cellVertices[0].y * (cellVertices[3].z - cellVertices[1].z) +
    cellVertices[1].y * (cellVertices[0].z - cellVertices[3].z) +
    cellVertices[3].y * (cellVertices[1].z - cellVertices[0].z)) +

    cellVertices[3].x * (cellVertices[0].y * (cellVertices[1].z - cellVertices[2].z) +
    cellVertices[1].y * (cellVertices[2].z - cellVertices[0].z) +
    cellVertices[2].y * (cellVertices[0].z - cellVertices[1].z));
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const //x -> ?
{
  Scalar invJacobian = Scalar(1.0) / GetCellDeformJacobian(cellVertices);
  Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices);
  Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices);
  Vector refZDerivatives = GetRefZDerivativesMulJacobian(cellVertices);
  
  return invJacobian * Vector(
    cellVertices[0].x * (cellVertices[3].y * cellVertices[2].z - cellVertices[2].y * cellVertices[3].z) +
    cellVertices[2].x * (cellVertices[0].y * cellVertices[3].z - cellVertices[3].y * cellVertices[0].z) +
    cellVertices[3].x * (cellVertices[2].y * cellVertices[0].z - cellVertices[0].y * cellVertices[2].z) +

    globalCoords.x * refXDerivatives.x + globalCoords.y * refXDerivatives.y + globalCoords.z * refXDerivatives.z,


    cellVertices[0].y * (cellVertices[3].x * cellVertices[1].z - cellVertices[1].x * cellVertices[3].z) +
    cellVertices[1].y * (cellVertices[0].x * cellVertices[3].z - cellVertices[3].x * cellVertices[0].z) +
    cellVertices[3].y * (cellVertices[1].x * cellVertices[0].z - cellVertices[0].x * cellVertices[1].z) +

    globalCoords.x * refYDerivatives.x + globalCoords.y * refYDerivatives.y + globalCoords.z * refYDerivatives.z,
    

    cellVertices[0].z * (cellVertices[2].x * cellVertices[1].y - cellVertices[1].x * cellVertices[2].y) +
    cellVertices[1].z * (cellVertices[0].x * cellVertices[2].y - cellVertices[2].x * cellVertices[0].y) +
    cellVertices[2].z * (cellVertices[1].x * cellVertices[0].y - cellVertices[0].x * cellVertices[1].y) +

    globalCoords.x * refZDerivatives.x + globalCoords.y * refZDerivatives.y + globalCoords.z * refZDerivatives.z
    );
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const // ? -> x
{
  return Vector(
    cellVertices[0].x + (cellVertices[1].x - cellVertices[0].x) * refCoords.x + 
                        (cellVertices[2].x - cellVertices[0].x) * refCoords.y + 
                        (cellVertices[3].x - cellVertices[0].x) * refCoords.z,
    cellVertices[0].y + (cellVertices[1].y - cellVertices[0].y) * refCoords.x + 
                        (cellVertices[2].y - cellVertices[0].y) * refCoords.y + 
                        (cellVertices[3].y - cellVertices[0].y) * refCoords.z,
    cellVertices[0].z + (cellVertices[1].z - cellVertices[0].z) * refCoords.x + 
                        (cellVertices[2].z - cellVertices[0].z) * refCoords.y + 
                        (cellVertices[3].z - cellVertices[0].z) * refCoords.z
    );
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GetRefXDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dξ/dx, dξ/dy, dξ/dz) * J
{
    return Vector(
      cellVertices[0].y * (cellVertices[2].z - cellVertices[3].z) + 
      cellVertices[2].y * (cellVertices[3].z - cellVertices[0].z) + 
      cellVertices[3].y * (cellVertices[0].z - cellVertices[2].z),
      cellVertices[0].x * (cellVertices[3].z - cellVertices[2].z) + 
      cellVertices[2].x * (cellVertices[0].z - cellVertices[3].z) + 
      cellVertices[3].x * (cellVertices[2].z - cellVertices[0].z),
      cellVertices[0].x * (cellVertices[2].y - cellVertices[3].y) + 
      cellVertices[2].x * (cellVertices[3].y - cellVertices[0].y) + 
      cellVertices[3].x * (cellVertices[0].y - cellVertices[2].y));
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GetRefYDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dη/dx, dη/dy, dη/dz) * J
{
    return Vector(
      cellVertices[0].y * (cellVertices[3].z - cellVertices[1].z) + 
      cellVertices[1].y * (cellVertices[0].z - cellVertices[3].z) + 
      cellVertices[3].y * (cellVertices[1].z - cellVertices[0].z),
      cellVertices[0].x * (cellVertices[1].z - cellVertices[3].z) + 
      cellVertices[1].x * (cellVertices[3].z - cellVertices[0].z) + 
      cellVertices[3].x * (cellVertices[0].z - cellVertices[1].z),
      cellVertices[0].x * (cellVertices[3].y - cellVertices[1].y) + 
      cellVertices[1].x * (cellVertices[0].y - cellVertices[3].y) + 
      cellVertices[3].x * (cellVertices[1].y - cellVertices[0].y));
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GetRefZDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dζ/dx, dζ/dy, dζ/dz) * J
{
    return Vector(
      cellVertices[0].y * (cellVertices[1].z - cellVertices[2].z) + 
      cellVertices[1].y * (cellVertices[2].z - cellVertices[0].z) + 
      cellVertices[2].y * (cellVertices[0].z - cellVertices[1].z),
      cellVertices[0].x * (cellVertices[2].z - cellVertices[1].z) + 
      cellVertices[1].x * (cellVertices[0].z - cellVertices[2].z) + 
      cellVertices[2].x * (cellVertices[1].z - cellVertices[0].z),
      cellVertices[0].x * (cellVertices[1].y - cellVertices[2].y) + 
      cellVertices[1].x * (cellVertices[2].y - cellVertices[0].y) + 
      cellVertices[2].x * (cellVertices[0].y - cellVertices[1].y));
}

template<typename FunctionSpace, typename System>
void  VolumeMesh<Space3, FunctionSpace, System>::GetRefDerivatives(
  Vector cellVertices[Space::NodesPerCell],
  Vector* refDerivatives) const
{
  Scalar invJacobian = Scalar(1.0) / fabs(GetCellDeformJacobian(cellVertices));
  refDerivatives[0] = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
  refDerivatives[1] = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;
  refDerivatives[2] = GetRefZDerivativesMulJacobian(cellVertices) * invJacobian;
}

template<typename FunctionSpace, typename System>
bool VolumeMesh<Space3, FunctionSpace, System>::IsCellRegular(IndexType cellIndex) const
{
  bool regularCell = true;
  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
  {
    IndexType interactionType = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;
    if (interactionType == IndexType(-1)) regularCell = false;
  }
  return regularCell;
}

template<typename FunctionSpace, typename System>
typename System::ValueType VolumeMesh<Space3, FunctionSpace, System>::GetFaceAverageSolution(IndexType cellIndex, IndexType faceNumber) const
{
  typename System::ValueType result(Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];

      result.values[valueIndex] += basisFunctionCoefficient * faceAverages[faceNumber].surfaceIntegral[functionIndex] * Scalar(2.0);
    }
  }

  return result;
}

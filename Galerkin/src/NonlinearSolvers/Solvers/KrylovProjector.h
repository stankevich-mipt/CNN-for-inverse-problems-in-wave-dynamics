#pragma once

#include "../NonlinearSystemSolver.h"
#include "../../Maths/MatrixMaths.h"

template <typename Scalar, typename IndexType>
struct KrylovProjectorSolver: public NonlinearSystemSolver<Scalar, IndexType>
{
  using NonlinearSystemSolver<Scalar, IndexType>::solverState;
  using NonlinearSystemSolver<Scalar, IndexType>::maxIterationsCount;
  using NonlinearSystemSolver<Scalar, IndexType>::tolerance;
  using NonlinearSystemSolver<Scalar, IndexType>::dimsCount;
  using NonlinearSystemSolver<Scalar, IndexType>::evalsCount;
  using NonlinearSystemSolver<Scalar, IndexType>::err;
  using NonlinearSystemSolver<Scalar, IndexType>::func;

  struct State
  {
    State(IndexType vectorsCount): 
      innerIteration(0), 
      innerState(Begin), 
      vectorsCount(vectorsCount)
    {
    }

    bool AdvanceStep()
    {
      bool phaseAdvanced = false;

      switch (innerState)
      {
        case Begin:
          innerState = FirstSyncRestore;
        break;
        case FirstSyncRestore:
          innerState = SecondSyncRestore;
        break;
        case SecondSyncRestore:
          innerIteration++;
          innerState = (innerIteration == vectorsCount) ? End : FirstSyncRestore;
        break;
        case End:
          innerIteration = 0;
          innerState = Begin;
          phaseAdvanced = true;
        break;
        default:
          assert(0);
        break;
      }
      return phaseAdvanced;
    }

    IndexType innerIteration;
    enum InnerState
    {
      Begin = 0,
      FirstSyncRestore,
      SecondSyncRestore,
      End
    } innerState;
    
    IndexType vectorsCount;
  };

  KrylovProjectorSolver(Scalar slopeMultiplier, IndexType iterationsCount, IndexType vectorsCount): state(vectorsCount)
  {
    this->slopeMultiplier = slopeMultiplier;
    this->vectorsCount    = vectorsCount;
  }

  virtual ~KrylovProjectorSolver()
  {
    delete [] currEstimation;
    delete [] nextEstimation;

    delete [] currValues;
    delete [] basisVectors;
    delete [] gramMatrix;
    delete [] gramTransposed;

    delete [] gramSymmetric;
    delete [] gramSymmetricInv;
    delete [] minimizerColumn;
    delete [] leastSquaresMult;

    delete [] basisCoords;

    delete [] directionalDerivatives;
    delete [] probeCoords;
  }

  virtual void SetFunction(NonlinearFunctionSystem<Scalar, IndexType> *func)
  {
    this->func = func;
    dimsCount = func->GetMaxDimentionsCount();

    currEstimation = new Scalar[dimsCount];
    nextEstimation = new Scalar[dimsCount];

    currValues = new Scalar[dimsCount];

    basisVectors = new Scalar[(vectorsCount + 1) * dimsCount];

    gramMatrix = new Scalar[(vectorsCount + 1) * vectorsCount];
    gramTransposed = new Scalar[vectorsCount * (vectorsCount + 1)];

    gramSymmetric = new Scalar[vectorsCount * vectorsCount];
    gramSymmetricInv = new Scalar[vectorsCount * vectorsCount];

    minimizerColumn = new Scalar[vectorsCount + 1];
    leastSquaresMult = new Scalar[vectorsCount * (vectorsCount + 1)];

    basisCoords = new Scalar[vectorsCount];

    directionalDerivatives = new Scalar[dimsCount];
    probeCoords = new Scalar[dimsCount];
  }

  void InitStep(const Scalar* initialGuess, Scalar* solution, const SolverState& solverState, IndexType maxIterationsCount, Scalar tolerance)
  {
    NonlinearSystemSolver<Scalar, IndexType>::InitStep(initialGuess, solution, solverState, maxIterationsCount, tolerance);
  }

  bool AdvancePhase(Scalar* solution, Scalar& err)
  { 
    const Scalar eps = Scalar(1e-4);

    switch (state.innerState)
    {
      case State::Begin:
      {
        evalsCount++;
        err = func->GetValue(solution, solverState, currValues);
      } break;

      case State::FirstSyncRestore: 
      {
        if (state.innerIteration == 0)
        {
          for(IndexType elementIndex = 0; elementIndex < vectorsCount * (vectorsCount + 1); elementIndex++)
          {
            gramMatrix[elementIndex] = 0;
          }

          len = 0;
          for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
          {
            len += sqr(currValues[coordIndex]);
          }

          len = sqrt(len);
          initialLen = len;

          for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
          {
            basisVectors[dimsCount * 0 + coordIndex] = -currValues[coordIndex] / len;
          }
        }
        evalsCount++;
        for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
        {
          probeCoords[coordIndex] = solution[coordIndex] + eps * basisVectors[dimsCount * state.innerIteration + coordIndex];
        }

        func->GetValue(probeCoords, solverState, directionalDerivatives);
      } break;

      case State::SecondSyncRestore:
      {
        for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
        {
          directionalDerivatives[coordIndex] = (directionalDerivatives[coordIndex] - currValues[coordIndex]) / eps;
        }

        for(IndexType vectorIndex = 0; vectorIndex < state.innerIteration + 1; vectorIndex++)
        {
          for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
          {
            gramMatrix[vectorIndex * vectorsCount + state.innerIteration] +=
              directionalDerivatives[coordIndex] * basisVectors[dimsCount * vectorIndex + coordIndex];
          }
        }

        //if(innerIteration + 1 >= vectorsCount) break;

        for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
        {
          basisVectors[(state.innerIteration + 1) * dimsCount + coordIndex] = directionalDerivatives[coordIndex];
        }
        for(IndexType vectorIndex = 0; vectorIndex < state.innerIteration + 1; vectorIndex++)
        {
          for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
          {
            basisVectors[(state.innerIteration + 1) * dimsCount + coordIndex] -=
              gramMatrix[vectorIndex * vectorsCount + state.innerIteration] * basisVectors[vectorIndex * dimsCount + coordIndex];
          }
        }

        len = 0;
        for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
        {
          len += sqr(basisVectors[(state.innerIteration + 1) * dimsCount + coordIndex]);
        }
        len = sqrt(len);
        for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
        {
          basisVectors[(state.innerIteration + 1) * dimsCount + coordIndex] /= len;
        }

        gramMatrix[(state.innerIteration + 1) * vectorsCount + state.innerIteration] = len;
      } break;

      case State::End:
      {
        for(IndexType vectorIndex = 0; vectorIndex < vectorsCount + 1; vectorIndex++)
        {
          minimizerColumn[vectorIndex] = 0;
        }
        minimizerColumn[0] = initialLen;

        MatrixTranspose(gramMatrix, gramTransposed, (vectorsCount + 1), vectorsCount);
        MatrixMulMatrix(gramTransposed, gramMatrix, gramSymmetric, vectorsCount, (vectorsCount + 1), vectorsCount);
        MatrixInverse(gramSymmetric, gramSymmetricInv, vectorsCount, Scalar(1e-5));
        MatrixMulMatrix(gramSymmetricInv, gramTransposed, leastSquaresMult, vectorsCount, vectorsCount, (vectorsCount + 1));
        MatrixMulMatrix<Scalar, IndexType>(leastSquaresMult, minimizerColumn, basisCoords, vectorsCount, vectorsCount + 1, 1);

        for(IndexType vectorIndex = 0; vectorIndex < vectorsCount; vectorIndex++)
        {
          for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
          {
            solution[coordIndex] += basisCoords[vectorIndex] * basisVectors[vectorIndex * dimsCount + coordIndex] * slopeMultiplier;
          }
        }
      } break;

      default:
        assert(0);
      break;
    }
    return state.AdvanceStep() && err < tolerance;
  }

  State state;
  IndexType vectorsCount;
  Scalar len;
  Scalar initialLen;

  Scalar slopeMultiplier;
  Scalar *currEstimation;
  Scalar *nextEstimation;

  Scalar *currValues;

  Scalar *basisVectors;
  Scalar *gramMatrix;
  Scalar *gramTransposed;

  Scalar *gramSymmetric;
  Scalar *gramSymmetricInv;

  Scalar *leastSquaresMult;

  Scalar *minimizerColumn;

  Scalar *basisCoords;
    
  Scalar *directionalDerivatives;
  Scalar *probeCoords;
};

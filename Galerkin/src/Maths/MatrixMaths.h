#pragma once
template<typename Scalar, typename IndexType>
void GaussSeidelIterate(const Scalar *matrix, const Scalar *rightVector, Scalar *currSolution, IndexType dimsCount)
{
  for(IndexType i = 0; i < dimsCount; i++)
  {
    Scalar val = rightVector[i];
    for(IndexType j = 0; j < dimsCount; j++)
    {
      if(i == j) continue;
      val -= matrix[i * dimsCount + j] * currSolution[j];
    }
    currSolution[i] = val / matrix[i * dimsCount + i];
  }
}

template<typename Scalar, typename IndexType, IndexType dimsCount>
void GaussSeidelSolve(const Scalar *matrix, const Scalar *rightVector, Scalar *solution, 
                      Scalar tolerance = Scalar(1e-5), IndexType maxIterations = 100)
{
  Scalar currSolution[dimsCount];
  for(IndexType i = 0; i < dimsCount; i++)
  {
    solution[i] = currSolution[i] = 0; 
  }
  Scalar err = tolerance * Scalar(2.0);
  for(IndexType iter = 0; ((iter < maxIterations) && (err > tolerance)); iter++)
  {
    GaussSeidelIterate<Scalar, IndexType>(matrix, rightVector, currSolution, dimsCount);
    err = 0;
    for(IndexType i = 0; i < dimsCount; i++)
    {
      err += fabs(currSolution[i] - solution[i]);
      solution[i] = currSolution[i];
    }
  }
}

template<typename Scalar, typename IndexType>
Scalar GetPivot(Scalar *matrix, Scalar *invMatrix, Scalar *rightVector, IndexType offset, const IndexType dimsCount)
{
  IndexType bestRow = offset;
  Scalar bestVal = 0;
  for(IndexType k = offset; k < dimsCount; k++)
  {
    Scalar val = fabs(matrix[k * dimsCount + offset]);
    if(val > bestVal)
    {
      bestVal = val;
      bestRow = k;
    }
  }

  if(bestRow > offset)
  {
    for(IndexType k = 0; k < dimsCount; k++)
    {
      Scalar t = matrix[offset * dimsCount + k];
      matrix[offset * dimsCount + k] = matrix[bestRow * dimsCount + k];
      matrix[bestRow * dimsCount + k] = t;
    }
    if(invMatrix)
    {
      for(IndexType k = 0; k < dimsCount; k++)
      {
        Scalar t = invMatrix[offset * dimsCount + k];
        invMatrix[offset * dimsCount + k] = invMatrix[bestRow * dimsCount + k];
        invMatrix[bestRow * dimsCount + k] = t;
      }
    }
    if(rightVector)
    {
      Scalar t = rightVector[bestRow];
      rightVector[bestRow] = rightVector[offset];
      rightVector[offset] = t;
    }
  }
  return matrix[offset * dimsCount + offset];
}        


template<typename Scalar, typename IndexType>
bool GaussSolveUnsafe(Scalar *matrix, Scalar *invMatrix, Scalar *rightVector, Scalar *solution, const IndexType dimsCount, const Scalar tolerance = Scalar(1e-5))
{
  //triangulation
  if(invMatrix)
  {
    for(IndexType i = 0; i < dimsCount; i++)
    {
      for(IndexType j = 0; j < dimsCount; j++)
      {
        invMatrix[i * dimsCount + j] = (i == j) ? Scalar(1.0) : Scalar(0.0);
      }
    }
  }

  for(IndexType offset = 0; offset < dimsCount; offset++)
  {
    //getting pivot
    IndexType bestRow = offset;
    Scalar bestVal = 0;
    for(IndexType k = offset; k < dimsCount; k++)
    {
      Scalar val = fabs(matrix[k * dimsCount + offset]);
      if(val > bestVal)
      {
        bestVal = val;
        bestRow = k;
      }
    }

    if(bestRow > offset)
    {
      for(IndexType k = 0; k < dimsCount; k++)
      {
        Scalar t = matrix[offset * dimsCount + k];
        matrix[offset * dimsCount + k] = matrix[bestRow * dimsCount + k];
        matrix[bestRow * dimsCount + k] = t;
      }
      if(invMatrix)
      {
        for(IndexType k = 0; k < dimsCount; k++)
        {
          Scalar t = invMatrix[offset * dimsCount + k];
          invMatrix[offset * dimsCount + k] = invMatrix[bestRow * dimsCount + k];
          invMatrix[bestRow * dimsCount + k] = t;
        }
      }
      if(rightVector)
      {
        Scalar t = rightVector[bestRow];
        rightVector[bestRow] = rightVector[offset];
        rightVector[offset] = t;
      }
    }
    //variable elimination
    Scalar diag = matrix[offset * dimsCount + offset];

    if(fabs(diag) < tolerance) return 0;

    for(IndexType j = 0; j < dimsCount; j++)
    {
      matrix[offset * dimsCount + j] /= diag;
    }

    if(invMatrix)
    {
      for(IndexType j = 0; j < dimsCount; j++)
      {
        invMatrix[offset * dimsCount + j] /= diag;
      }
    }

    if(rightVector)
    {
      rightVector[offset] /= diag;
    }

    for(IndexType j = offset + 1; j < dimsCount; j++)
    {
      Scalar mult = matrix[j * dimsCount + offset];
      for(IndexType k = 0; k < dimsCount; k++)
      {
        matrix[j * dimsCount + k] -= mult * matrix[offset * dimsCount + k];
      }
      if(rightVector)
      {
        rightVector[j] -= mult * rightVector[offset];
      }
      if(invMatrix)
      {
        for(IndexType k = 0; k < dimsCount; k++)
        {
          invMatrix[j * dimsCount + k] -= mult * invMatrix[offset * dimsCount + k];
        }
      }
    }
  }
  //solution substitution
  if(rightVector && solution)
  {
    for(IndexType i = 0; i < dimsCount; i++)
    {
      IndexType row = dimsCount - i - 1;
      Scalar sum = 0;
      for(IndexType j = row + 1; j < dimsCount; j++)
      {
        sum += matrix[row * dimsCount + j] * solution[j];
      }
      solution[row] = (rightVector[row] - sum) / matrix[row * dimsCount + row];
    }
  }

  for(IndexType i = 0; i < dimsCount; i++)
  {
    IndexType column = dimsCount - i - 1;
    IndexType row = column;
    assert(matrix[row * dimsCount + column] != 0);
    Scalar mult = Scalar(1.0) / matrix[row * dimsCount + column];

    for(IndexType j = 0; j < dimsCount; j++)
    {
      matrix[row * dimsCount + j] *= mult;
    }

    if(invMatrix)
    {
      for(IndexType j = 0; j < dimsCount; j++)
      {
        invMatrix[row * dimsCount + j] *= mult; 
      }
    }
    if(rightVector)
    {
      rightVector[row] *= mult;
    }

    for(IndexType eliminatingRow = 0; eliminatingRow < row; eliminatingRow++)
    {
      Scalar mult = matrix[eliminatingRow * dimsCount + column];
      for(IndexType j = 0; j < dimsCount; j++)
      {
        matrix[eliminatingRow * dimsCount + j] -= matrix[row * dimsCount + j] * mult;
      }

      if(invMatrix)
      {
        for(IndexType j = 0; j < dimsCount; j++)
        {
          invMatrix[eliminatingRow * dimsCount + j] -= invMatrix[row * dimsCount + j] * mult;
        }
      }

      if(rightVector)
      {
        rightVector[eliminatingRow] -= rightVector[row] * mult;
      }
    }
  }
  return 1;
}

/*template<typename Scalar, typename IndexType, IndexType dimsCount>
void GaussSolve(Scalar *matrix, Scalar *invMatrix, Scalar *rightVector, Scalar *solution, const Scalar tolerance = Scalar(1e-5))
{
  Scalar tmpMatrix[dimsCount * dimsCount];
  Scalar tmpVector[dimsCount];
  for(IndexType i = 0; i < dimsCount; i++)
  {
    for(IndexType j = 0; j < dimsCount; j++)
    {
      tmpMatrix[i * dimsCount + j] = matrix[i * dimsCount + j];
    }
    tmpVector[i] = rightVector[i];
  }
  GaussSolveUnsafe<Scalar, IndexType>(tmpMatrix, tmpVector, solution, dimsCount, tolerance);
}*/

template<typename Scalar, typename IndexType>
void GaussSolve(Scalar *matrix, Scalar *rightVector, Scalar *solution, const IndexType dimsCount, const Scalar tolerance = Scalar(1e-5))
{
  Scalar *tmpMatrix = new Scalar[dimsCount * dimsCount];
  Scalar *tmpVector = new Scalar[dimsCount];
  for(IndexType i = 0; i < dimsCount; i++)
  {
    for(IndexType j = 0; j < dimsCount; j++)
    {
      tmpMatrix[i * dimsCount + j] = matrix[i * dimsCount + j];
    }
    if(rightVector)
    {
      tmpVector[i] = rightVector[i];
    }
  }
  if(!rightVector)
  {
    delete[] tmpVector;
    tmpVector = 0;
  }

  GaussSolveUnsafe<Scalar, IndexType>(tmpMatrix, 0, tmpVector, solution, dimsCount, tolerance);
  delete [] tmpMatrix;
  if(tmpVector)
  {
    delete [] tmpVector;
  }
}

template<typename Scalar, typename IndexType>
void MatrixCopy(Scalar *srcMatrix, Scalar *dstMatrix, IndexType height, IndexType width)
{
  std::copy(srcMatrix, srcMatrix + height * width, dstMatrix);
  /* for(IndexType i = 0; i < height; i++)
  {
    for(IndexType j = 0; j < width; j++)
    {
      dstMatrix[width * i + j] = srcMatrix[width * i + j];
    }
  } */
}

template<typename Scalar, typename IndexType>
bool MatrixInverse(Scalar *matrix, Scalar *invMatrix, const IndexType dimsCount, const Scalar tolerance = Scalar(1e-9))
{
  std::vector<Scalar> tmpMatrix(dimsCount * dimsCount);
  MatrixCopy(matrix, tmpMatrix.data(), dimsCount, dimsCount);
  return GaussSolveUnsafe<Scalar, IndexType>(tmpMatrix.data(), invMatrix, 0, 0, dimsCount, tolerance);
}

template<typename Scalar, typename IndexType>
void MatrixTranspose(Scalar *srcMatrix, Scalar *dstMatrix, IndexType height, IndexType width)
{
  for(IndexType i = 0; i < height; i++)
  {
    for(IndexType j = 0; j < width; j++)
    {
      dstMatrix[height * j + i] = srcMatrix[width * i + j];
    }
  }
}

template<typename Scalar, typename IndexType>
void MatrixMulMatrix(Scalar *matrix0, Scalar *matrix1, Scalar *resMatrix, IndexType matrix0Height, IndexType matrix0Width, IndexType matrix1Width)
{
  memset(resMatrix, 0, matrix0Height * matrix1Width * sizeof(Scalar));

  register Scalar t;
  int i, j, k;
  for (i = 0; i < matrix0Height; ++i)
  {
    for (k = 0; k < matrix0Width; ++k)
    {
      t = matrix0[matrix0Width * i + k];
      for (j = 0; j < matrix1Width; ++j)
      {
        resMatrix[i * matrix1Width + j] += t * matrix1[k * matrix1Width + j];
      }
    }
  }
}

template<typename Scalar, typename IndexType>
void MatrixMulTransponsedMatrix(Scalar *matrix0, Scalar *matrix1, Scalar *resMatrix, IndexType matrix0Height, IndexType matrix0Width, IndexType matrix1Width)
{
  register Scalar t;
  int i, j, k;
  for (i = 0; i < matrix0Height; ++i)
  {
    for (j = 0; j < matrix1Width; ++j)
    {
      t = 0;
      for (k = 0; k < matrix0Width; ++k)
      {
        t += matrix0[matrix0Width * i + k] * matrix1[j * matrix1Width + k];
      }
      resMatrix[i * matrix1Width + j] = t;
    }
  }
}

template<typename Scalar, typename IndexType>
void MatrixMulVector(Scalar *matrix, Scalar *vec, Scalar *resVec, IndexType matrixHeight, IndexType matrixWidth)
{
  for(IndexType i = 0; i < matrixHeight; i++)
  {
    Scalar sum = 0;
    for(IndexType j = 0; j < matrixWidth; j++)
    {
      sum += matrix[matrixWidth * i + j] * vec[j];
    }

    resVec[i] = sum;
  }
}

template<typename Scalar, typename IndexType>
void MatrixAddMatrix(Scalar *matrix0, Scalar *matrix1, IndexType height, IndexType width)
{
  for(IndexType i = 0; i < height; i++)
  {
    for(IndexType j = 0; j < width; j++)
    {
      matrix0[width * i + j] += matrix1[width * i + j];
    }
  }
}

template<typename Scalar, typename IndexType>
void MatrixAddMatrix(Scalar *matrix0, Scalar *matrix1, Scalar *resMatrix, IndexType height, IndexType width)
{
  for(IndexType i = 0; i < height; i++)
  {
    for(IndexType j = 0; j < width; j++)
    {
      resMatrix[width * i + j] = matrix0[width * i + j] + matrix1[width * i + j];
    }
  }
}

template<typename Scalar, typename IndexType>
void MatrixSubMatrix(Scalar *matrix0, Scalar *matrix1, IndexType height, IndexType width)
{
  for(IndexType i = 0; i < height; i++)
  {
    for(IndexType j = 0; j < width; j++)
    {
      matrix0[width * i + j] -= matrix1[width * i + j];
    }
  }
}

template<typename Scalar, typename IndexType>
void MatrixSubMatrix(Scalar *matrix0, Scalar *matrix1, Scalar *resMatrix, IndexType height, IndexType width)
{
  for(IndexType i = 0; i < height; i++)
  {
    for(IndexType j = 0; j < width; j++)
    {
      resMatrix[width * i + j] = matrix0[width * i + j] - matrix1[width * i + j];
    }
  }
}

template<typename Scalar, typename IndexType>
void MatrixMulScalar(Scalar *matrix, Scalar mult, IndexType height, IndexType width)
{
  for(IndexType i = 0; i < height; i++)
  {
    for(IndexType j = 0; j < width; j++)
    {
      matrix[width * i + j] *= mult;
    }
  }
}

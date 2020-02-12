#pragma once
#include "../../../Maths/MatrixMaths.h"
#include "../../../Maths/Spaces.h"
#include "BaseSpace.h"

template<typename Space, int Order>
struct PolynomialSpace;

template <int Order>
struct PolynomialSpace<Space2, Order>: public BaseSpace<Space2, Order>
{
  SPACE2_TYPEDEFS

  using BaseSpace<Space2, Order>::order;
  using BaseSpace<Space2, Order>::functionsCount;

  Scalar functionMult;

  PolynomialSpace()
  {
    functionMult = 1.0;
    ComputeFunctionIndexByPows();
    ComputeCoordPowers();
    ComputeMaxFunctionValues();
    ComputeVolumeIntegrals();
  }

  Scalar GetBasisFunctionValue(Vector2 point, IndexType functionIndex) const
  {
    /*{
      switch(functionIndex)
      {
      case 0: return Scalar(1.0); break;
      case 1: return point.x; break;
      case 2: return point.y; break;
      case 3: return point.z; break;
      }
    }*/
    return GetPolynomeValue(GetCoordPowers(functionIndex), point) * functionMult;
  }

  template<typename Function, int ValuesCount>
  void Decompose(Function func, Scalar *coords)
  {
    Scalar correlations[functionsCount * ValuesCount]; // correlations == inner products
    for(IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
    {
      for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
      {
        correlations[valueIndex * functionsCount + functionIndex] = 0;
        coords      [valueIndex * functionsCount + functionIndex] = 0;
      }
    }
    
    Scalar step = Scalar(1e-2);
    for(Scalar x = step / Scalar(2.0); x <= Scalar(1.0); x += step)
    {
      for(Scalar y = step / Scalar(2.0); y <= Scalar(1.0); y += step)
      {
        if(x + y <= Scalar(1.0))
        {
          Scalar values[ValuesCount];
          func(Vector2(x, y), values);
          for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
          {
            Scalar basisFunctionValue = GetBasisFunctionValue(Vector2(x, y), functionIndex);
            for(IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
            {
              correlations[valueIndex * functionsCount + functionIndex] += basisFunctionValue * values[valueIndex] * step * step;
            }
          }
        }
      }
    }

    for(IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
    {
      for(IndexType i = 0; i < functionsCount; i++)
      {
        for(IndexType j = 0; j < functionsCount; j++)
        {
          coords[valueIndex * functionsCount + i] += correlations[valueIndex * functionsCount + j] * cellVolumeIntegralsInv[j * functionsCount + i];
        }
      }
    }
  }

  Polynomial<Scalar, IndexType, 2> GetBasisPolynomial(IndexType functionIndex) const
  {
    IndexVector polynomialPows = GetCoordPowers(functionIndex);

    typename Polynomial<Scalar, IndexType, 2>::Pows pows;
    pows.pows[0] = polynomialPows.x;
    pows.pows[1] = polynomialPows.y;

    Polynomial<Scalar, IndexType, 2> res;
    res.AddTerm(pows, Scalar(1.0));
    return res;
  }

  Scalar GetBasisFunctionMaxValue(IndexType functionIndex) const
  {
    /* float maxValue = 0;
    float testVal;
    testVal = GetBasisFunctionValue(Vector2(0, 0), functionIndex); */
    return 0;
  }
  Scalar GetMaxFunctionValue(IndexType functionIndex)
  {
    return maxFunctionValues[functionIndex];
  }

  Scalar GetBasisFunctionDerivative(Vector2 point, IndexVector derivatives, IndexType functionIndex)
  {
    /*{
      switch(functionIndex)
      {
      case 0: return ((derivatives.x > 0) || (derivatives.y > 0) || (derivatives.z > 0)) ? 0 : Scalar(1.0); break;
      case 1:
        {
          if(derivatives.x == 1)
          {
            if(derivatives.y == 0 && derivatives.z == 0)
              return Scalar(1.0);
            else
              return 0;
          }else
            return 0;
        } break;
      case 2: 
        {
          if(derivatives.y == 1)
          {
            if(derivatives.x == 0 && derivatives.z == 0)
              return Scalar(1.0);
            else
              return 0;
          }else
            return 0;
        } break;
      case 3:
        {
          if(derivatives.z == 1)
          {
            if(derivatives.x == 0 && derivatives.y == 0)
              return Scalar(1.0);
            else
              return 0;
          }else
            return 0;
        } break;
      }
    }
    return 0;*/

    IndexVector resPows;
    IndexType coef;
    GetPolynomeDerivative(GetCoordPowers(functionIndex), derivatives, resPows, coef);
    return coef * GetPolynomeValue(resPows, point) * functionMult;
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex)
  {
    IndexVector pows = GetCoordPowers(functionIndex);

    Scalar analyticalRes = ComputePolynomialIntegral(pows) * functionMult;
    return analyticalRes;
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ф1, Ф2>
  {
    /*Scalar res = 0;

    Scalar step = Scalar(1e-3);
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        if(x + y < Scalar(1.0))
        {
          res += GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                 GetBasisFunctionValue(Vector2(x, y), functionIndex1) *
                 step * step;
        }
      }
    }*/
    //return res;

    IndexVector pows0 = GetCoordPowers(functionIndex0);
    IndexVector pows1 = GetCoordPowers(functionIndex1);

    IndexVector pows = pows0 + pows1;
    Scalar analyticalRes = ComputePolynomialIntegral(pows) * functionMult * functionMult;
    /*if(fabs(analyticalRes - res) > 1e-3)
    {
      printf("\nErr: %f %f i0: %d i1: %d", analyticalRes, res, functionIndex0, functionIndex1);
    }*/
    return analyticalRes;
  }

  Vector2 ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<dФ1/dx, Ф2>, <dФ1/dy, Ф2>, <dФ1/dz, Ф2>}
  {
    Vector2 res(0, 0);
    IndexVector pows0 = GetCoordPowers(functionIndex0);
    IndexVector pows1 = GetCoordPowers(functionIndex1);

    if(pows1.x > 0)
    {
      res.x = Scalar(pows1.x) * ComputePolynomialIntegral((pows0 + pows1 + IndexVector(-1, 0)));
    }
    if(pows1.y > 0)
    {
      res.y = Scalar(pows1.y) * ComputePolynomialIntegral((pows0 + pows1 + IndexVector(0, -1)));
    }
    return res * functionMult * functionMult;

    /*Vector2 res(0, 0);

    Scalar step = Scalar(1e-3);
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        if(x + y < Scalar(1.0))
        {
          res.x +=  GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                    GetBasisFunctionDerivative(Vector2(x, y), IndexVector(1, 0), functionIndex1) *
                    step * step;
          res.y +=  GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                    GetBasisFunctionDerivative(Vector2(x, y), IndexVector(0, 1), functionIndex1) *
                    step * step;
        }
      }
    }
    return res;*/
  }


public:
  Scalar Pow(Scalar a, Scalar b) const
  {
    if(fabs(b) < Scalar(1e-5)) return Scalar(1.0);
    if(fabs(a) < Scalar(1e-5)) return 0;
    return exp(b * log(a));
  }

  void GetPolynomeDerivative(IndexVector pows, IndexVector derivatives, IndexVector &resPows, IndexType &coef) const
  {
    coef = 1;
    if(derivatives.x > pows.x)
    {
      coef = 0;
      resPows = IndexVector(0, 0);
      return;
    }else
    {
      resPows.x = pows.x;
      for(IndexType i = 0; i < derivatives.x; i++)
      {
        coef *= resPows.x;
        resPows.x--;
      }
    }

    if(derivatives.y > pows.y)
    {
      coef = 0;
      resPows = IndexVector(0, 0);
      return;
    }else
    {
      resPows.y = pows.y;
      for(IndexType i = 0; i < derivatives.y; i++)
      {
        coef *= resPows.y;
        resPows.y--;
      }
    }
  }

  Scalar GetPolynomeValue(IndexVector pows, Vector2 point) const
  {
    return Pow(point.x, Scalar(pows.x)) * Pow(point.y, Scalar(pows.y));
  }


  void ComputeMaxFunctionValues()
  {
    for(int functionIndex = 0; functionIndex < functionsCount; functionIndex++)
    {
      Scalar maxValue = 0;
      Scalar testValue;

      testValue = fabs(GetBasisFunctionValue(Vector2(0, 0), functionIndex));
      if(testValue > maxValue) maxValue = testValue;

      testValue = fabs(GetBasisFunctionValue(Vector2(1, 0), functionIndex));
      if(testValue > maxValue) maxValue = testValue;

      testValue = fabs(GetBasisFunctionValue(Vector2(0, 1), functionIndex));
      if(testValue > maxValue) maxValue = testValue;

      maxFunctionValues[functionIndex] = maxValue;
    }
  }

  void ComputeCoordPowers()
  {
    IndexType num = 0;
    for(int xPow = 0; xPow <= Order; xPow++)
    {
      for(int yPow = 0; yPow <= Order; yPow++)
      {
        if(xPow + yPow <= Order)
        {
          num++;
          functionPows[num - 1] = IndexVector(xPow, yPow);
        }
      }
    }
  }

  IndexVector GetCoordPowers(IndexType functionIndex) const
  {
    return functionPows[functionIndex];
  }

  void ComputeFunctionIndexByPows()
  {
    IndexType num = 0;
    for(int xPow = 0; xPow <= Order; xPow++)
    {
      for(int yPow = 0; yPow <= Order; yPow++)
      {
        functionIndices[xPow][yPow] = -1;
        if(xPow + yPow <= Order)
        {
          num++;
          functionIndices[xPow][yPow] = num - 1;
        }
      }
    }
  }

  IndexType GetFunctionIndexByPows(IndexVector pows)
  {
    return functionIndices[pows.x][pows.y];
  }

  static Scalar Factorial(IndexType num)
  {
    Scalar res(1.0);
    for(IndexType i = 1; i <= num; i++)
    {
      res *= Scalar(i);
    }
    return res;
  }

  Scalar ComputePolynomialIntegral(IndexVector pows) const
  {
    return 
      Scalar(Factorial(pows.x) * Factorial(pows.y)) /
      Scalar(Factorial(pows.x + pows.y + 2));
  }

  void ComputeVolumeIntegrals()
  {
    for(IndexType functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
    {
      for(IndexType functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
      {
        cellVolumeIntegrals[functionIndex0 * functionsCount + functionIndex1] =
          ComputeCellVolumeIntegral       (functionIndex1, functionIndex0);
      }
    }
    MatrixInverse<Scalar, IndexType>(cellVolumeIntegrals, cellVolumeIntegralsInv, functionsCount);
  }

  IndexType functionIndices[Order + 1][Order + 1];
  IndexVector functionPows[functionsCount];
  Scalar   maxFunctionValues[functionsCount];
  Scalar   cellVolumeIntegrals    [functionsCount * functionsCount];
  Scalar   cellVolumeIntegralsInv [functionsCount * functionsCount];
};

template<int Order>
struct PolynomialSpace<Space3, Order>: public BaseSpace<Space3, Order>
{
  SPACE3_TYPEDEFS

  using BaseSpace<Space3, Order>::order;
  using BaseSpace<Space3, Order>::functionsCount;

  PolynomialSpace()
  {
    ComputeFunctionIndexByPows();
    ComputeCoordPowers();
    ComputeVolumeIntegrals();
  }

  Scalar GetBasisFunctionValue(Vector3 point, IndexType functionIndex) const
  {
    return GetPolynomeValue(GetCoordPowers(functionIndex), point);
  }

  template<typename Function, int ValuesCount>
  void Decompose(Function func, Scalar *coords)
  {
    Scalar correlations[functionsCount * ValuesCount];
    for (IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
    {
      for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
      {
        correlations[valueIndex * functionsCount + functionIndex] = 0;
        coords[valueIndex * functionsCount + functionIndex] = 0;
      }
    }

    Scalar step = Scalar(1e-2);
    for (Scalar x = step / Scalar(2.0); x <= Scalar(1.0); x += step)
    {
      for (Scalar y = step / Scalar(2.0); y <= Scalar(1.0); y += step)
      {
        for (Scalar z = step / Scalar(2.0); z <= Scalar(1.0); z += step)
        {
          if (x + y + z <= Scalar(1.0))
          {
            Scalar values[ValuesCount];
            func(Vector3(x, y, z), values);

            for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
            {
              Scalar basisFunctionValue = GetBasisFunctionValue(Vector3(x, y, z), functionIndex);

              for (IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
              {
                correlations[valueIndex * functionsCount + functionIndex] += basisFunctionValue * values[valueIndex] * step * step * step;
              }
            }
          }
        }
      }
    }

    for (IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
    {
      for (IndexType i = 0; i < functionsCount; i++)
      {
        for (IndexType j = 0; j < functionsCount; j++)
        {
          coords[valueIndex * functionsCount + i] += correlations[valueIndex * functionsCount + j] * cellVolumeIntegralsInv[j * functionsCount + i];
        }
      }
    }
  }

  Polynomial<Scalar, IndexType, 3> GetBasisPolynomial(IndexType functionIndex) const
  {
    IndexVector polynomialPows = GetCoordPowers(functionIndex);

    typename Polynomial<Scalar, IndexType, 3>::Pows pows;
    pows.pows[0] = polynomialPows.x;
    pows.pows[1] = polynomialPows.y;
    pows.pows[2] = polynomialPows.z;

    Polynomial<Scalar, IndexType, 3> res;
    res.AddTerm(pows, Scalar(1.0));
    return res;
  }

  Scalar GetBasisFunctionDerivative(Vector point, IndexVector derivatives, IndexType functionIndex)
  {
    IndexVector resPows;
    IndexType coef;
    GetPolynomeDerivative(GetCoordPowers(functionIndex), derivatives, resPows, coef);
    return coef * GetPolynomeValue(resPows, point);
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ô1, Ô2>
  {
    IndexVector pows0 = GetCoordPowers(functionIndex0);
    IndexVector pows1 = GetCoordPowers(functionIndex1);

    IndexVector pows = pows0 + pows1;
    return ComputePolynomialIntegral(pows);
    /*

    Scalar res = 0;

    Scalar step = 5e-3f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    for(Scalar z = 0; z < Scalar(1.0); z += step)
    {
    if(x + y + z < Scalar(1.0))
    {
    res += GetBasisFunctionValue(Vector3(x, y, z), functionIndex0) *
    GetBasisFunctionValue(Vector3(x, y, z), functionIndex1) *
    step * step * step;
    }
    }
    }
    }
    return res;*/
  }

  Vector3 ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<Ô2, dÔ1/dx>, <Ô2, dÔ1/dy>, <Ô2, dÔ1/dz>}
  {
    Vector3 res(0, 0, 0);
    IndexVector pows0 = GetCoordPowers(functionIndex0);
    IndexVector pows1 = GetCoordPowers(functionIndex1);

    if (pows1.x > 0)
    {
      res.x = Scalar(pows1.x) * ComputePolynomialIntegral((pows0 + pows1 + IndexVector(-1, 0, 0)));
    }
    if (pows1.y > 0)
    {
      res.y = Scalar(pows1.y) * ComputePolynomialIntegral((pows0 + pows1 + IndexVector(0, -1, 0)));
    }
    if (pows1.z > 0)
    {
      res.z = Scalar(pows1.z) * ComputePolynomialIntegral((pows0 + pows1 + IndexVector(0, 0, -1)));
    }
    return res;
    /*
    Vector3 res(0, 0, 0);

    Scalar step = 1e-2f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    for(Scalar z = 0; z < Scalar(1.0); z += step)
    {
    if(x + y + z < Scalar(1.0))
    {
    res.x +=  GetBasisFunctionValue(Vector3(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector3(x, y, z), IndexVector(1, 0, 0), functionIndex1) *
    step * step * step;
    res.y +=  GetBasisFunctionValue(Vector3(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector3(x, y, z), IndexVector(0, 1, 0), functionIndex1) *
    step * step * step;
    res.z +=  GetBasisFunctionValue(Vector3(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector3(x, y, z), IndexVector(0, 0, 1), functionIndex1) *
    step * step * step;
    }
    }
    }
    }
    return res;*/
  }


private:
  Scalar Pow(Scalar a, Scalar b) const
  {
    if (fabs(b) < Scalar(1e-5)) return Scalar(1.0);
    if (fabs(a) < Scalar(1e-5)) return 0;
    return exp(b * log(a));
  }

  void GetPolynomeDerivative(IndexVector pows, IndexVector derivatives, IndexVector &resPows, IndexType &coef) const
  {
    coef = 1;
    if (derivatives.x > pows.x)
    {
      coef = 0;
      resPows = IndexVector(0, 0, 0);
      return;
    } else
    {
      resPows.x = pows.x;
      for (IndexType i = 0; i < derivatives.x; i++)
      {
        coef *= resPows.x;
        resPows.x--;
      }
    }

    if (derivatives.y > pows.y)
    {
      coef = 0;
      resPows = IndexVector(0, 0, 0);
      return;
    } else
    {
      resPows.y = pows.y;
      for (IndexType i = 0; i < derivatives.y; i++)
      {
        coef *= resPows.y;
        resPows.y--;
      }
    }

    if (derivatives.z > pows.z)
    {
      coef = 0;
      resPows = IndexVector(0, 0, 0);
      return;
    } else
    {
      resPows.z = pows.z;
      for (IndexType i = 0; i < derivatives.z; i++)
      {
        coef *= resPows.z;
        resPows.z--;
      }
    }
  }

  Scalar GetPolynomeValue(IndexVector pows, Vector3 point) const
  {
    return Pow(point.x, Scalar(pows.x)) * Pow(point.y, Scalar(pows.y)) * Pow(point.z, Scalar(pows.z));
  }


  void ComputeCoordPowers()
  {
    IndexType num = 0;
    for (int xPow = 0; xPow <= Order; xPow++)
    {
      for (int yPow = 0; yPow <= Order; yPow++)
      {
        for (int zPow = 0; zPow <= Order; zPow++)
        {
          if (xPow + yPow + zPow <= Order)
          {
            num++;
            functionPows[num - 1] = IndexVector(xPow, yPow, zPow);
          }
        }
      }
    }
  }

  IndexVector GetCoordPowers(IndexType functionIndex) const
  {
    return functionPows[functionIndex];
  }

  void ComputeFunctionIndexByPows()
  {
    IndexType num = 0;
    for (int xPow = 0; xPow <= Order; xPow++)
    {
      for (int yPow = 0; yPow <= Order; yPow++)
      {
        for (int zPow = 0; zPow <= Order; zPow++)
        {
          functionIndices[xPow][yPow][zPow] = -1;
          if (xPow + yPow + zPow <= Order)
          {
            num++;
            functionIndices[xPow][yPow][zPow] = num - 1;
          }
        }
      }
    }
  }

  IndexType GetFunctionIndexByPows(IndexVector pows)
  {
    return functionIndices[pows.x][pows.y][pows.z];
  }

  static Scalar Factorial(IndexType num)
  {
    Scalar res(1.0);
    for (IndexType i = 1; i <= num; i++)
    {
      res *= Scalar(i);
    }
    return res;
  }

  Scalar ComputePolynomialIntegral(IndexVector pows) const
  {
    return
      Scalar(Factorial(pows.x) * Factorial(pows.y) * Factorial(pows.z)) /
      Scalar(Factorial(pows.x + pows.y + pows.z + 3));
  }

  void ComputeVolumeIntegrals()
  {
    for (IndexType functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
    {
      for (IndexType functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
      {
        cellVolumeIntegrals[functionIndex0 * functionsCount + functionIndex1] =
          ComputeCellVolumeIntegral(functionIndex1, functionIndex0);
      }
    }
    MatrixInverse<Scalar, IndexType>(cellVolumeIntegrals, cellVolumeIntegralsInv, functionsCount);
  }

  IndexType functionIndices[Order + 1][Order + 1][Order + 1];
  IndexVector functionPows[functionsCount];
  Scalar   cellVolumeIntegrals[functionsCount * functionsCount];
  Scalar   cellVolumeIntegralsInv[functionsCount * functionsCount];
};

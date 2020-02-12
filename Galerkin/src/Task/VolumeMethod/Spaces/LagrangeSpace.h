#pragma once
#include "../../../Maths/Spaces.h"
#include "BaseSpace.h"

template<typename Space, int Order>
struct LagrangeSpace;

template<int Order>
struct LagrangeSpace<Space2, Order>: public BaseSpace<Space2, Order>
{
  SPACE2_TYPEDEFS

  using BaseSpace<Space2, Order>::order;
  using BaseSpace<Space2, Order>::functionsCount;

  Scalar functionMult;

  typedef Polynomial<Scalar, IndexType, 2> PolynomialType;

  LagrangeSpace(): BaseSpace<Space2, Order>()
  {
    ComputeCoordPowers();
  }

  Polynomial<Scalar, IndexType, 2> GetBasisPolynomial(IndexType functionIndex) const
  {
    return GetBasisFunctionPolynomial(functionIndex);
    /*Vector2i polynomialPows = GetCoordPowers(functionIndex);

    typename Polynomial<Scalar, IndexType, 2>::Pows pows;
    pows.pows[0] = polynomialPows.x;
    pows.pows[1] = polynomialPows.y;

    Polynomial<Scalar, IndexType, 2> res;
    res.AddTerm(pows, Scalar(1.0));
    return res;*/
  }


  PolynomialType GetBasisFunctionPolynomial(IndexType functionIndex) const
  {
    Vector localCoordGradients[3];

    Vector vertices[3];
    vertices[0] = Vector(0, 0);
    vertices[1] = Vector(Scalar(1.0), 0);
    vertices[2] = Vector(0, Scalar(1.0));

    Scalar square = Scalar(0.5) * (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]);
    Scalar invSquare = Scalar(1.0) / square;

    Vector dirs[3];
    dirs[0] = vertices[2] - vertices[1];
    dirs[1] = vertices[0] - vertices[2];
    dirs[2] = vertices[1] - vertices[0];

    localCoordGradients[0] = Vector(-dirs[0].y, dirs[0].x) * Scalar(0.5) * invSquare;
    localCoordGradients[1] = Vector(-dirs[1].y, dirs[1].x) * Scalar(0.5) * invSquare;
    localCoordGradients[2] = Vector(-dirs[2].y, dirs[2].x) * Scalar(0.5) * invSquare;

    PolynomialType localCoordPolynomials[3];

    typename PolynomialType::Pows defPows;
    defPows.pows[0] = 0;
    defPows.pows[1] = 0;

    localCoordPolynomials[2].AddTerm(defPows, Scalar(1.0));
    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      typename PolynomialType::Pows termPows;

      termPows.pows[0] = termPows.pows[1] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, -localCoordGradients[coordIndex] * vertices[coordIndex + 1]);

      termPows.pows[0] = 1;
      termPows.pows[1] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].x);

      termPows.pows[0] = 0;
      termPows.pows[1] = 1;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].y);

      localCoordPolynomials[2] -= localCoordPolynomials[coordIndex];
    }

    IndexType coordPowers[3];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar denom(1.0);

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        denom /= Scalar(termNumber + 1) / Scalar(Order);
      }
    }

    PolynomialType resPolynomial;

    resPolynomial.AddTerm(defPows, denom);

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        PolynomialType localTerm;

        typename PolynomialType::Pows localPows;
        localPows.pows[0] = 0;
        localPows.pows[1] = 0;
        localTerm.AddTerm(localPows, -Scalar(termNumber) / Scalar(Order));

        localTerm += localCoordPolynomials[coordIndex];

        resPolynomial = resPolynomial * localTerm;
      }
    }

    return resPolynomial;
  }

  Vector GetBasisPoint(IndexType functionIndex) const
  {
    IndexType coordPowers[3];
    GetCoordPowers(functionIndex, coordPowers);

    Vector vertices[3];
    vertices[0] = Vector(0, 0);
    vertices[1] = Vector(Scalar(1.0), 0);
    vertices[2] = Vector(0, Scalar(1.0));

    return vertices[0] 
      + (vertices[1] - vertices[0]) * Scalar(coordPowers[1]) / Scalar(Order)
      + (vertices[2] - vertices[0]) * Scalar(coordPowers[2]) / Scalar(Order);
  }

  Scalar GetBasisFunctionValue(const Vector& point, IndexType functionIndex) const
  {
    IndexType coordPowers[3];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar localCoords[3];
    GetLocalCoords(point, localCoords);

    return GetPolynomeValue(coordPowers, localCoords);

    /*PolynomialType polynomial = GetBasisFunctionPolynomial(functionIndex);
    Scalar polynomialRes = polynomial.GetValue(&(point.x));

    if(fabs(normalRes - polynomialRes) > 1e-5) printf("poo");
    return polynomialRes;*/
  }

  Scalar GetBasisFunctionMaxValue(IndexType functionIndex)
  {
    return Scalar(1.0);
  }

  

public:
  void GetLocalCoords(const Vector& point, Scalar* localCoords) const
  {
    Vector vertices[3];
    vertices[0] = Vector(0, 0);
    vertices[1] = Vector(Scalar(1.0), 0);
    vertices[2] = Vector(0, Scalar(1.0));
    Scalar square = Scalar(0.5) * (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]);
    Scalar invSquare = Scalar(1.0) / square;

    localCoords[0] = ((vertices[2] - vertices[1]) ^ (point - vertices[1])) * invSquare * Scalar(0.5);
    localCoords[1] = ((vertices[0] - vertices[2]) ^ (point - vertices[2])) * invSquare * Scalar(0.5);
    localCoords[2] = ((vertices[1] - vertices[0]) ^ (point - vertices[0])) * invSquare * Scalar(0.5);
  }


  Scalar GetPolynomeValue(IndexType* pows, Scalar* localCoords) const
  {
    Scalar res(1.0);

    for(IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res /= Scalar(power + 1) / Scalar(Order);
      }
    }

    // Scalar denom = res;

    for(IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res *= localCoords[coordIndex] - Scalar(power) / Scalar(Order);
      }
    }

    // Scalar baseRes = res;
    return res;
  }

  void ComputeCoordPowers()
  {
    IndexType functionIndex = 0;
    for(IndexType i = 0; i <= Order; i++)
    {
      for(IndexType j = 0; j <= Order - i; j++)
      {
        IndexType k = Order - i - j;
        functionPows[functionIndex][0] = i;
        functionPows[functionIndex][1] = j;
        functionPows[functionIndex][2] = k;
        functionIndex++;
      }
    }
  }

  void GetCoordPowers(IndexType functionIndex, IndexType* pows) const
  {
    for(IndexType i = 0; i < 3; i++)
    {
      pows[i] = functionPows[functionIndex][i];
    }
  }
  IndexType functionPows[functionsCount][3];
};

template<int Order>
struct LagrangeSpace<Space3, Order>: public BaseSpace<Space3, Order>
{
  SPACE3_TYPEDEFS

  using BaseSpace<Space3, Order>::order;
  using BaseSpace<Space3, Order>::functionsCount;

  Scalar functionMult;

  typedef Polynomial<Scalar, IndexType, 3> PolynomialType;

  LagrangeSpace()
  {
    ComputeCoordPowers();
  }

  Polynomial<Scalar, IndexType, 3> GetBasisPolynomial(IndexType functionIndex) const
  {
    return GetBasisFunctionPolynomial(functionIndex);
  }


  PolynomialType GetBasisFunctionPolynomial(IndexType functionIndex) const
  {

    Vector localCoordGradients[4];

    Vector vertices[4];
    vertices[0] = Vector(0, 0, 0);
    vertices[1] = Vector(Scalar(1), 0, 0);
    vertices[2] = Vector(0, Scalar(1), 0);
    vertices[3] = Vector(0, 0, Scalar(1));

    Scalar volume = MixedProduct(vertices[1] - vertices[0], vertices[2] - vertices[0], vertices[3] - vertices[0]) / Scalar(6.0);
    // Scalar invVolume = Scalar(1.0) / volume;


    localCoordGradients[0] = (vertices[3] - vertices[1]) ^ (vertices[2] - vertices[1]) / (Scalar(6.0) * volume);
    localCoordGradients[1] = (vertices[2] - vertices[0]) ^ (vertices[3] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[2] = (vertices[3] - vertices[0]) ^ (vertices[1] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[3] = (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]) / (Scalar(6.0) * volume);

    PolynomialType localCoordPolynomials[4];

    PolynomialType::Pows defPows;
    defPows.pows[0] = 0;
    defPows.pows[1] = 0;
    defPows.pows[2] = 0;

    localCoordPolynomials[3].AddTerm(defPows, Scalar(1.0));
    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      PolynomialType::Pows termPows;

      termPows.pows[0] = termPows.pows[1] = termPows.pows[2] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, -localCoordGradients[coordIndex] * vertices[coordIndex + 1]);

      termPows.pows[0] = 1;
      termPows.pows[1] = 0;
      termPows.pows[2] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].x);

      termPows.pows[0] = 0;
      termPows.pows[1] = 1;
      termPows.pows[2] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].y);

      termPows.pows[0] = 0;
      termPows.pows[1] = 0;
      termPows.pows[2] = 1;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].z);

      localCoordPolynomials[3] -= localCoordPolynomials[coordIndex];
    }

    IndexType coordPowers[4];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar denom(1.0);

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        denom /= Scalar(termNumber + 1) / Scalar(Order);
      }
    }

    PolynomialType resPolynomial;

    resPolynomial.AddTerm(defPows, denom);

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        PolynomialType localTerm;

        PolynomialType::Pows localPows;
        localPows.pows[0] = 0;
        localPows.pows[1] = 0;
        localPows.pows[2] = 0;
        localTerm.AddTerm(localPows, -Scalar(termNumber) / Scalar(Order));

        localTerm += localCoordPolynomials[coordIndex];

        resPolynomial = resPolynomial * localTerm;
      }
    }

    return resPolynomial;
  }

  Vector GetBasisPoint(IndexType functionIndex) const
  {
    IndexType coordPowers[4];
    GetCoordPowers(functionIndex, coordPowers);

    Vector vertices[4];
    vertices[0] = Vector(0, 0, 0);
    vertices[1] = Vector(Scalar(1), 0, 0);
    vertices[2] = Vector(0, Scalar(1), 0);
    vertices[3] = Vector(0, 0, Scalar(1));

    return vertices[0]
      + (vertices[1] - vertices[0]) * Scalar(coordPowers[1]) / Scalar(Order)
      + (vertices[2] - vertices[0]) * Scalar(coordPowers[2]) / Scalar(Order)
      + (vertices[3] - vertices[0]) * Scalar(coordPowers[3]) / Scalar(Order);

  }

  Scalar GetBasisFunctionValue(const Vector& point, IndexType functionIndex) const
  {
    IndexType coordPowers[4];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar localCoords[4];
    GetLocalCoords(point, localCoords);

    return GetPolynomeValue(coordPowers, localCoords);

    /*PolynomialType polynomial = GetBasisFunctionPolynomial(functionIndex);
    Scalar polynomialRes = polynomial.GetValue(&(point.x));

    if(fabs(normalRes - polynomialRes) > 1e-5) printf("poo");
    return polynomialRes;*/
  }

  Scalar GetBasisFunctionMaxValue(IndexType functionIndex)
  {
    return Scalar(1.0);
  }

public:

  void GetLocalCoords(const Vector& point, Scalar *localCoords) const
  {
    Vector vertices[4];
    vertices[0] = Vector(0, 0, 0);
    vertices[1] = Vector(Scalar(1), 0, 0);
    vertices[2] = Vector(0, Scalar(1), 0);
    vertices[3] = Vector(0, 0, Scalar(1));

    Scalar volume = MixedProduct(vertices[1] - vertices[0], vertices[2] - vertices[0], vertices[3] - vertices[0]) / Scalar(6.0);
    // Scalar invVolume = Scalar(1.0) / volume;

    Vector localCoordGradients[4];

    localCoordGradients[0] = (vertices[3] - vertices[1]) ^ (vertices[2] - vertices[1]) / (Scalar(6.0) * volume);
    localCoordGradients[1] = (vertices[2] - vertices[0]) ^ (vertices[3] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[2] = (vertices[3] - vertices[0]) ^ (vertices[1] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[3] = (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]) / (Scalar(6.0) * volume);

    localCoords[0] = localCoordGradients[0] * (point - vertices[1]);
    localCoords[1] = localCoordGradients[1] * (point - vertices[2]);
    localCoords[2] = localCoordGradients[2] * (point - vertices[3]);
    localCoords[3] = localCoordGradients[3] * (point - vertices[0]);
  }

  Scalar GetPolynomeValue(IndexType *pows, Scalar *localCoords) const
  {
    Scalar res(1.0);

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res /= Scalar(power + 1) / Scalar(Order);
      }
    }

    // Scalar denom = res;

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res *= localCoords[coordIndex] - Scalar(power) / Scalar(Order);
      }
    }

    // Scalar baseRes = res;

    return res;
  }

  void ComputeCoordPowers()
  {
    IndexType functionIndex = 0;
    for (IndexType i = 0; i <= Order; i++)
    {
      for (IndexType j = 0; j <= Order - i; j++)
      {
        for (IndexType k = 0; k <= Order - i - j; k++)
        {
          IndexType l = Order - i - j - k;
          functionPows[functionIndex][0] = i;
          functionPows[functionIndex][1] = j;
          functionPows[functionIndex][2] = k;
          functionPows[functionIndex][3] = l;
          functionIndex++;
        }
      }
    }
  }

  void GetCoordPowers(IndexType functionIndex, IndexType *pows) const
  {
    for (IndexType i = 0; i < 4; i++)
    {
      pows[i] = functionPows[functionIndex][i];
    }
  }
  IndexType functionPows[functionsCount][4];
};

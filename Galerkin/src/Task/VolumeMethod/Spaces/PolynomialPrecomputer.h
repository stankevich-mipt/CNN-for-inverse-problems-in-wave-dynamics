#pragma once

#include "../Cell.h"
#include "Decomposer.h" 

template<typename Space, typename Decomposer>
struct PolynomialPrecomputer;

template<typename Decomposer>
struct PolynomialPrecomputer<Space2, Decomposer>: public Decomposer
{
  SPACE2_TYPEDEFS
  typedef Space2 Space;

  using Decomposer::functionsCount;
  using Decomposer::order;

  PolynomialPrecomputer()
  {
  }

  Scalar GetBasisFunctionDerivative(Vector point, IndexVector derivatives, IndexType functionIndex)
  {
    Scalar coords[2] = {point.x, point.y};
    return GetBasisFunctionDerivative(derivatives, functionIndex).GetValue(coords);
  }

  Polynomial<Scalar, IndexType, 2> GetBasisFunctionDerivative(IndexVector derivatives, IndexType functionIndex)
  {
    IndexType derivativeDegrees[2] = {derivatives.x, derivatives.y};
    Polynomial<Scalar, IndexType, 2> function = this->GetBasisPolynomial(functionIndex);
    return function.Differentiate(derivativeDegrees);
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex) 
  {
    Polynomial<Scalar, IndexType, 2> function = this->GetBasisPolynomial(functionIndex);
    return function.ComputeSubspaceIntegral(2);
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ф1, Ф2>
  {
    Polynomial<Scalar, IndexType, 2> function0 = this->GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 2> function1 = this->GetBasisPolynomial(functionIndex1);

    Polynomial<Scalar, IndexType, 2> correlation = function0 * function1;
    return correlation.ComputeSubspaceIntegral(2);
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
            res += GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                   GetBasisFunctionValue(Vector(x, y, z), functionIndex1) *
                   step * step * step;
          }
        }
      }
    }
    return res;*/
  }

  Vector ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<Ф2, dФ1/dx>, <Ф2, dФ1/dy>, <Ф2, dФ1/dz>}
  {
    Polynomial<Scalar, IndexType, 2> function0 = this->GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 2> function1 = this->GetBasisPolynomial(functionIndex1);

    IndexType derivatives[2];

    derivatives[0] = 1;
    derivatives[1] = 0;
    Polynomial<Scalar, IndexType, 2> dfunction1dx = function1.Differentiate(derivatives);

    derivatives[0] = 0;
    derivatives[1] = 1;
    Polynomial<Scalar, IndexType, 2> dfunction1dy = function1.Differentiate(derivatives);

    return Vector((function0 * dfunction1dx).ComputeSubspaceIntegral(2), (function0 * dfunction1dy).ComputeSubspaceIntegral(2));
    /*
    Vector res(0, 0, 0);

    Scalar step = 1e-2f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        for(Scalar z = 0; z < Scalar(1.0); z += step)
        {
          if(x + y + z < Scalar(1.0))
          {
            res.x +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                      GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(1, 0, 0), functionIndex1) *
                      step * step * step;
            res.y +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                      GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 1, 0), functionIndex1) *
                      step * step * step;
            res.z +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                      GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 0, 1), functionIndex1) *
                      step * step * step;
          }
        }
      }
    }
    return res;*/
  }

  Scalar ComputeOutgoingFlux(IndexType refEdgeNumber, IndexType basisFunction0, IndexType basisFunction1)
  {
    /*Scalar step = Scalar(1e-5) * 10;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      res += 
        GetBasisFunctionValue(TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber, x), basisFunction1) *
        GetBasisFunctionValue(TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber, x), basisFunction0) * step;
    } 
    return res;*/

    Polynomial<Scalar, IndexType, 2> function0 = this->GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 2> function1 = this->GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 1> edgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(refEdgeNumber, edgeToRef);

    typename Polynomial<Scalar, IndexType, 1>::Pows pows;

    Polynomial<Scalar, IndexType, 1> intToEdge[1];
    pows.pows[0] = 1;
    intToEdge[0].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 1> intToRef[2];

    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      intToRef[coordIndex] = edgeToRef[coordIndex].Substitute(intToEdge);
    }

    Polynomial<Scalar, IndexType, 1> intToSrcFunc;
    Polynomial<Scalar, IndexType, 1> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToRef);
    intToDstFunc = function0.Substitute(intToRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(1);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }

  Scalar ComputeIncomingFlux(
    IndexType refEdgeNumber, IndexType correspondingEdgeNumber,
    IndexType basisFunction0, IndexType basisFunction1)
  {
    /*if(refEdgeNumber == 0 && correspondingEdgeNumber == 1 && basisFunction0 == 0 && basisFunction1 == 1)
    {
      int pp = 1;
    }
    Scalar step = Scalar(1e-5) * 10;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      Scalar srcEdgeCoords = x;
      Scalar dstEdgeCoords = TriCell::GetCorrespondingEdgeCoord<GeomSpace>(srcEdgeCoords);

      Vector2 srcRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber,           srcEdgeCoords);
      Vector2 dstRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(correspondingEdgeNumber, dstEdgeCoords);

      res += 
        GetBasisFunctionValue(dstRefCoords, basisFunction0) *
        GetBasisFunctionValue(srcRefCoords, basisFunction1) * step;
    } */
    //return res;

    Polynomial<Scalar, IndexType, 2> function0 = this->GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 2> function1 = this->GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 1> srcEdgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(refEdgeNumber          , srcEdgeToRef);

    Polynomial<Scalar, IndexType, 1> dstEdgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(correspondingEdgeNumber, dstEdgeToRef);


    typename Polynomial<Scalar, IndexType, 1>::Pows pows;

    Polynomial<Scalar, IndexType, 1> intToSrcEdge[1];
    pows.pows[0] = 1;
    intToSrcEdge[0].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 1> srcEdgeToDstEdge[1];
    Cell<Space2>::GetCorrespondingEdgeTransitionPolynomials(srcEdgeToDstEdge);

    Polynomial<Scalar, IndexType, 1> intToDstEdge[1];
    intToDstEdge[0] = srcEdgeToDstEdge[0].Substitute(intToSrcEdge);

    Polynomial<Scalar, IndexType, 1> intToSrcRef[2];
    Polynomial<Scalar, IndexType, 1> intToDstRef[2];

    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      intToSrcRef[coordIndex] = srcEdgeToRef[coordIndex].Substitute(intToSrcEdge);
      intToDstRef[coordIndex] = dstEdgeToRef[coordIndex].Substitute(intToDstEdge);
    }

    Polynomial<Scalar, IndexType, 1> intToSrcFunc;
    Polynomial<Scalar, IndexType, 1> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToSrcRef);
    intToDstFunc = function0.Substitute(intToDstRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(1);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }

  Scalar ComputeEdgeFlux(IndexType refEdgeNumber, IndexType basisFunction)
  {
    Polynomial<Scalar, IndexType, 2> function = this->GetBasisPolynomial(basisFunction);

    Polynomial<Scalar, IndexType, 1> edgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(refEdgeNumber, edgeToRef);

    typename Polynomial<Scalar, IndexType, 1>::Pows pows;

    Polynomial<Scalar, IndexType, 1> intToEdge[1];
    intToEdge[0] = Polynomial<Scalar, IndexType, 1>("x");

    Polynomial<Scalar, IndexType, 1> intToRef[2];

    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      intToRef[coordIndex] = edgeToRef[coordIndex].Substitute(intToEdge);
    }

    Polynomial<Scalar, IndexType, 1> intToFunc;

    intToFunc = function.Substitute(intToRef);

    Scalar integrationResult = intToFunc.ComputeSubspaceIntegral(1);
    return integrationResult;
  }
};

template<typename Decomposer>
struct PolynomialPrecomputer<Space3, Decomposer>: public Decomposer
{
  SPACE3_TYPEDEFS
  typedef Space3 Space;

  using Decomposer::functionsCount;
  using Decomposer::order;

  PolynomialPrecomputer()
  {
  }

  Scalar GetBasisFunctionDerivative(Vector point, IndexVector derivatives, IndexType functionIndex)
  {
    Scalar coords[3] = {point.x , point.y , point.z};
    return GetBasisFunctionDerivative(derivatives, functionIndex).GetValue(coords);
  }

  Polynomial<Scalar, IndexType, 3> GetBasisFunctionDerivative(IndexVector derivatives, IndexType functionIndex)
  {
    Polynomial<Scalar, IndexType, 3> function = this->GetBasisPolynomial(functionIndex);
    IndexType derivativeDegrees[3] = {derivatives.x , derivatives.y , derivatives.z};
    return function.Differentiate(derivativeDegrees);
  }


  Scalar ComputeCellVolumeIntegral(IndexType functionIndex)
  {
    Polynomial<Scalar, IndexType, 3> function = this->GetBasisPolynomial(functionIndex);
    return function.ComputeSubspaceIntegral(3);
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ô1, Ô2>
  {
    Polynomial<Scalar, IndexType, 3> function0 = this->GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 3> function1 = this->GetBasisPolynomial(functionIndex1);

    Polynomial<Scalar, IndexType, 3> correlation = function0 * function1;
    return correlation.ComputeSubspaceIntegral(3);

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
    res += GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionValue(Vector(x, y, z), functionIndex1) *
    step * step * step;
    }
    }
    }
    }
    return res;*/
  }

  Vector ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<Ô2, dÔ1/dx>, <Ô2, dÔ1/dy>, <Ô2, dÔ1/dz>}
  {
    Polynomial<Scalar, IndexType, 3> function0 = this->GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 3> function1 = this->GetBasisPolynomial(functionIndex1);

    IndexType derivatives[3];

    derivatives[0] = 1;
    derivatives[1] = 0;
    derivatives[2] = 0;
    Polynomial<Scalar, IndexType, 3> dfunction1dx = function1.Differentiate(derivatives);

    derivatives[0] = 0;
    derivatives[1] = 1;
    derivatives[2] = 0;
    Polynomial<Scalar, IndexType, 3> dfunction1dy = function1.Differentiate(derivatives);

    derivatives[0] = 0;
    derivatives[1] = 0;
    derivatives[2] = 1;
    Polynomial<Scalar, IndexType, 3> dfunction1dz = function1.Differentiate(derivatives);

    return Vector((function0 * dfunction1dx).ComputeSubspaceIntegral(3),
      (function0 * dfunction1dy).ComputeSubspaceIntegral(3),
      (function0 * dfunction1dz).ComputeSubspaceIntegral(3));
    /*
    Vector res(0, 0, 0);

    Scalar step = 1e-2f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    for(Scalar z = 0; z < Scalar(1.0); z += step)
    {
    if(x + y + z < Scalar(1.0))
    {
    res.x +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(1, 0, 0), functionIndex1) *
    step * step * step;
    res.y +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 1, 0), functionIndex1) *
    step * step * step;
    res.z +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 0, 1), functionIndex1) *
    step * step * step;
    }
    }
    }
    }
    return res;*/
  }

  Scalar ComputeOutgoingFlux(IndexType refFaceNumber, IndexType basisFunction0, IndexType basisFunction1)
  {
    /*Scalar step = 1e-2f;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    if(x + y < Scalar(1.0))
    {
    res +=
    functionSpace.GetBasisFunctionValue(TetraCell::GetSurfaceRefCoords(refFaceNumber, Vector2(x, y)), basisFunction1) *
    functionSpace.GetBasisFunctionValue(TetraCell::GetSurfaceRefCoords(refFaceNumber, Vector2(x, y)), basisFunction0) * step * step;
    }
    }
    }
    return res;*/

    Polynomial<Scalar, IndexType, 3> function0 = this->GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 3> function1 = this->GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 2> faceToRef[3];
    Cell<Space3>::GetFaceRefTransitionPolynomials(refFaceNumber, faceToRef);

    Polynomial<Scalar, IndexType, 2>::Pows pows;

    Polynomial<Scalar, IndexType, 2> intToFace[2];

    pows.pows[0] = 1;
    pows.pows[1] = 0;
    intToFace[0].AddTerm(pows, 1);

    pows.pows[0] = 0;
    pows.pows[1] = 1;
    intToFace[1].AddTerm(pows, 1);
    
    Polynomial<Scalar, IndexType, 2> intToRef[3];

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      intToRef[coordIndex] = faceToRef[coordIndex].Substitute(intToFace);
    }

    Polynomial<Scalar, IndexType, 2> intToSrcFunc;
    Polynomial<Scalar, IndexType, 2> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToRef);
    intToDstFunc = function0.Substitute(intToRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(2);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }

  Scalar ComputeIncomingFlux(
    IndexType refFaceNumber, IndexType correspondingFaceNumber, IndexType orientation,
    IndexType basisFunction0, IndexType basisFunction1)
  {
    /*Scalar step = 1e-2f;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    if(x + y < Scalar(1.0))
    {
    Vector2 srcSurfaceCoords = Vector2(x, y);
    Vector2 dstSurfaceCoords = TetraCell::GetCorrespondingSurfaceCoords(srcSurfaceCoords, orientation);

    Vector srcRefCoords   = TetraCell::GetSurfaceRefCoords(refFaceNumber,           srcSurfaceCoords);
    Vector dstRefCoords   = TetraCell::GetSurfaceRefCoords(correspondingFaceNumber, dstSurfaceCoords);

    res +=
    functionSpace.GetBasisFunctionValue(dstRefCoords, basisFunction0) *
    functionSpace.GetBasisFunctionValue(srcRefCoords, basisFunction1) * step * step;
    }
    }
    }
    return res;*/
    /*if(refEdgeNumber == 0 && correspondingEdgeNumber == 1 && basisFunction0 == 0 && basisFunction1 == 1)
    {
    int pp = 1;
    }
    Scalar step = Scalar(1e-5) * 10;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    Scalar srcEdgeCoords = x;
    Scalar dstEdgeCoords = TriCell::GetCorrespondingEdgeCoord<GeomSpace>(srcEdgeCoords);

    Vector2 srcRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber,           srcEdgeCoords);
    Vector2 dstRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(correspondingEdgeNumber, dstEdgeCoords);

    res +=
    GetBasisFunctionValue(dstRefCoords, basisFunction0) *
    GetBasisFunctionValue(srcRefCoords, basisFunction1) * step;
    } */
    //return res;

    Polynomial<Scalar, IndexType, 3> function0 = this->GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 3> function1 = this->GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 2> srcFaceToRef[3];
    Cell<Space3>::GetFaceRefTransitionPolynomials(refFaceNumber, srcFaceToRef);

    Polynomial<Scalar, IndexType, 2> dstFaceToRef[3];
    Cell<Space3>::GetFaceRefTransitionPolynomials(correspondingFaceNumber, dstFaceToRef);


    Polynomial<Scalar, IndexType, 2>::Pows pows;

    Polynomial<Scalar, IndexType, 2> intToSrcFace[2];

    pows.pows[0] = 1;
    pows.pows[1] = 0;
    intToSrcFace[0].AddTerm(pows, 1);

    pows.pows[0] = 0;
    pows.pows[1] = 1;
    intToSrcFace[1].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 2> srcFaceToDstFace[2];
    Cell<Space3>::GetCorrespondingFaceTransitionPolynomials(orientation, srcFaceToDstFace);

    Polynomial<Scalar, IndexType, 2> intToDstFace[2];
    intToDstFace[0] = srcFaceToDstFace[0].Substitute(intToSrcFace);
    intToDstFace[1] = srcFaceToDstFace[1].Substitute(intToSrcFace);

    Polynomial<Scalar, IndexType, 2> intToSrcRef[3];
    Polynomial<Scalar, IndexType, 2> intToDstRef[3];

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      intToSrcRef[coordIndex] = srcFaceToRef[coordIndex].Substitute(intToSrcFace);
      intToDstRef[coordIndex] = dstFaceToRef[coordIndex].Substitute(intToDstFace);
    }

    Polynomial<Scalar, IndexType, 2> intToSrcFunc;
    Polynomial<Scalar, IndexType, 2> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToSrcRef);
    intToDstFunc = function0.Substitute(intToDstRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(2);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }
  Scalar ComputeFaceFlux(IndexType refFaceNumber, IndexType basisFunction)
  {
    Polynomial<Scalar, IndexType, 3> function = this->GetBasisPolynomial(basisFunction);

    Polynomial<Scalar, IndexType, 2> faceToRef[3];
    Cell<Space3>::GetFaceRefTransitionPolynomials(refFaceNumber, faceToRef);

    Polynomial<Scalar, IndexType, 2> intToFace[2];
    intToFace[0] = Polynomial<Scalar, IndexType, 2>("x");
    intToFace[1] = Polynomial<Scalar, IndexType, 2>("y");

    Polynomial<Scalar, IndexType, 2> intToRef[3];

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      intToRef[coordIndex] = faceToRef[coordIndex].Substitute(intToFace);
    }

    Polynomial<Scalar, IndexType, 2> intToFunc;

    intToFunc = function.Substitute(intToRef);

    Scalar integrationResult = intToFunc.ComputeSubspaceIntegral(2);
    return integrationResult;
  }
};


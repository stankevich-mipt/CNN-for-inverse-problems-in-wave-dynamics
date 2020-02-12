#pragma once

#include "../../../src/Maths/Spaces.h"
#include <math.h>

template <typename Space>
struct Cell;

template <>
struct Cell<Space2>
{
  SPACE2_TYPEDEFS

  static Vector GetNode(IndexType nodeNumber)
  {
    switch (nodeNumber)
    {
      case 0: return Vector(0, 0); break;
      case 1: return Vector(1, 0); break;
      case 2: return Vector(0, 1); break;
      default: assert(0); return Vector(-1, -1); break;
    }
  }

  static Vector GetEdgeRefCoords(IndexType refEdgeNumber, Scalar edgeCoord)
  {
    switch(refEdgeNumber)
    {
      case 0:
      {
        return Vector(edgeCoord, 0);
      } break;
      case 1:
      {
        return Vector(Scalar(1.0) - edgeCoord, edgeCoord);
      } break;
      case 2:
      {
        return Vector(0, Scalar(1.0) - edgeCoord);
      } break;
    }
    assert(0);
    return Vector(-1, -1);
  }

  static Scalar GetEdgeRefLen(IndexType refEdgeNumber)
  {
    return (GetEdgeRefCoords(refEdgeNumber, 1) - GetEdgeRefCoords(refEdgeNumber, 0)).Len();
  }

  static void GetEdgeRefTransitionPolynomials(IndexType refEdgeNumber, Polynomial<Scalar, IndexType, 1>* coords)
  {
    typedef Polynomial<Scalar, IndexType, 1>::Pows PowsType;
    PowsType pows;

    switch(refEdgeNumber)
    {
      case 0:
      {
        pows.pows[0] = 1;
        coords[0].AddTerm(pows, 1);

        pows.pows[0] = 0;
        coords[1].AddTerm(pows, 0);
        return;
        //return Vector(edgeCoord, 0);
      } break;
      case 1:
      {
        pows.pows[0] = 0;
        coords[0].AddTerm(pows, 1);
        pows.pows[0] = 1;
        coords[0].AddTerm(pows, -1);

        pows.pows[0] = 1;
        coords[1].AddTerm(pows, 1);
        return;
        //return Vector(Scalar(1.0) - edgeCoord, edgeCoord);
      } break;
      case 2:
      {
        pows.pows[0] = 0;
        coords[0].AddTerm(pows, 0);

        pows.pows[0] = 0;
        coords[1].AddTerm(pows, 1);
        pows.pows[0] = 1;
        coords[1].AddTerm(pows, -1);
        return;
        //return Vector(0, Scalar(1.0) - edgeCoord);
      } break;
    }
    assert(0);
    return;
  }

  static Scalar GetCorrespondingEdgeCoord(Scalar srcEdgeCoord)
  {
    return Scalar(1.0) - srcEdgeCoord;
  }

  static void GetCorrespondingEdgeTransitionPolynomials(/*IndexType orientation,*/
    Polynomial<Scalar, IndexType, 1>* dstCoords)
  {
    typedef Polynomial<Scalar, IndexType, 1>::Pows PowsType;

    PowsType pows;

    pows.pows[0] = 0;
    dstCoords[0].AddTerm(pows, 1);

    pows.pows[0] = 1;
    dstCoords[0].AddTerm(pows, -1);

    /*switch(refEdgeNumber)
    {
      case 0:
      {
        pows[0] = 1;
        coords[0].AddTerm(pows, 1);

        pows[0] = 0;
        coords[1].AddTerm(pows, 0);
        return;

        //return Vector(edgeCoord, 0);
      }
      break;
      case 1:
      {
        pows[0] = 0;
        coords[0].AddTerm(pows, 1);
        pows[0] = 1;
        coords[0].AddTerm(pows, -1);

        pows[0] = 1;
        coords[1].AddTerm(pows, 1);
        return;
  //      return Vector(Scalar(1.0) - edgeCoord, edgeCoord);
      }
      break;
      case 2:
      {
        pows[0] = 0;
        coords[0].AddTerm(pows, 0);

        pows[0] = 0;
        coords[1].AddTerm(pows, 1);
        pows[0] = 1;
        coords[1].AddTerm(pows, -1);
        return;
  //      return Vector(0, Scalar(1.0) - edgeCoord);
      }
      break;
    }
    assert(0);
    return;*/
  }

  static Scalar GetRefCellVolume()
  {
    return Scalar(0.5);
  }
};

template <>
struct Cell<Space3>
{
  SPACE3_TYPEDEFS

  static Vector GetNode(IndexType nodeNumber)
  {
    switch (nodeNumber)
    {
      case 0: return Vector(0, 0, 0); break;
      case 1: return Vector(1, 0, 0); break;
      case 2: return Vector(0, 1, 0); break;
      case 3: return Vector(0, 0, 1); break;
      default: assert(0); return Vector(-1, -1, -1); break;
    }
  }

  static Vector GetSurfaceRefCoords(IndexType refFaceNumber, Vector2 surfaceCoords)
  {
    switch (refFaceNumber)
    {
      case 0:
      {
        return Vector(surfaceCoords.y, surfaceCoords.x, 0);
      } break;
      case 1:
      {
        return Vector(surfaceCoords.x, 0, surfaceCoords.y);
      } break;
      case 2:
      {
        return Vector(0, surfaceCoords.y, surfaceCoords.x);
      } break;
      case 3:
      {
        return Vector(Scalar(1.0) - surfaceCoords.x - surfaceCoords.y, surfaceCoords.x, surfaceCoords.y);
      } break;
    }
    assert(0);
    return Vector(-1, -1, -1);
  }

  static void GetFaceRefTransitionPolynomials(IndexType refFaceNumber, Polynomial<Scalar, IndexType, 2> *coords)
  {
    typedef Polynomial<Scalar, IndexType, 2>::Pows PowsType;
    PowsType pows;

    switch (refFaceNumber)
    {
    case 0:
    {
      pows.pows[0] = 0;
      pows.pows[1] = 1;
      coords[0].AddTerm(pows, 1);

      pows.pows[0] = 1;
      pows.pows[1] = 0;
      coords[1].AddTerm(pows, 1);

      pows.pows[0] = 0;
      pows.pows[1] = 0;
      coords[2].AddTerm(pows, 0);
      return;

      //return Vector2(edgeCoord, 0);
      //return typename Space::Vector3(surfaceCoords.y, surfaceCoords.x, 0);
    }
    break;
    case 1:
    {
      pows.pows[0] = 1;
      pows.pows[1] = 0;
      coords[0].AddTerm(pows, 1);

      pows.pows[0] = 0;
      pows.pows[1] = 0;
      coords[1].AddTerm(pows, 0);

      pows.pows[0] = 0;
      pows.pows[1] = 1;
      coords[2].AddTerm(pows, 1);
      return;
      //      return Vector2(Scalar(1.0) - edgeCoord, edgeCoord);
      //      return typename Space::Vector3(surfaceCoords.x, 0, surfaceCoords.y);
    }
    break;
    case 2:
    {
      pows.pows[0] = 0;
      pows.pows[1] = 0;
      coords[0].AddTerm(pows, 0);

      pows.pows[0] = 0;
      pows.pows[1] = 1;
      coords[1].AddTerm(pows, 1);

      pows.pows[0] = 1;
      pows.pows[1] = 0;
      coords[2].AddTerm(pows, 1);
      return;
      //      return Vector2(0, Scalar(1.0) - edgeCoord);
      //return typename Space::Vector3(0, surfaceCoords.y, surfaceCoords.x);
    }
    break;
    case 3:
    {
      pows.pows[0] = 0;
      pows.pows[1] = 0;
      coords[0].AddTerm(pows, 1);
      pows.pows[0] = 1;
      pows.pows[1] = 0;
      coords[0].AddTerm(pows, -1);
      pows.pows[0] = 0;
      pows.pows[1] = 1;
      coords[0].AddTerm(pows, -1);

      pows.pows[0] = 1;
      pows.pows[1] = 0;
      coords[1].AddTerm(pows, 1);

      pows.pows[0] = 0;
      pows.pows[1] = 1;
      coords[2].AddTerm(pows, 1);
      return;
      //      return Vector2(0, Scalar(1.0) - edgeCoord);
      /*pows.pows[0] = 0;
      coords[0].AddTerm(pows, 1);
      pows.pows[0] = 1;
      coords[0].AddTerm(pows, -1);

      pows.pows[0] = 1;
      coords[1].AddTerm(pows, 1);
      return;
      return typename Space::Vector3(Scalar(1.0) - surfaceCoords.x - surfaceCoords.y, surfaceCoords.x, surfaceCoords.y);*/
    }
    break;
    }
    assert(0);
    return;
  }

  static Vector2 GetCorrespondingSurfaceCoords(Vector2 srcSurfaceCoords, IndexType orientation)
  {
    switch (orientation)
    {
      case 0:
      {
        return Vector2(srcSurfaceCoords.y, srcSurfaceCoords.x);
      } break;
      case 1:
      {
        return Vector2(Scalar(1.0) - srcSurfaceCoords.x - srcSurfaceCoords.y, srcSurfaceCoords.y);
      } break;
      case 2:
      {
        return Vector2(srcSurfaceCoords.x, Scalar(1.0) - srcSurfaceCoords.x - srcSurfaceCoords.y);
      } break;
    }
    assert(0);
    return Vector2(Scalar(-1), Scalar(-1));
  }

  static void GetCorrespondingFaceTransitionPolynomials(IndexType orientation,
    Polynomial<Scalar, IndexType, 2> *dstCoords)
  {
    typedef Polynomial<Scalar, IndexType, 2>::Pows PowsType;
    PowsType pows;

    switch (orientation)
    {
    case 0:
    {
      pows.pows[0] = 0;
      pows.pows[1] = 1;
      dstCoords[0].AddTerm(pows, 1);

      pows.pows[0] = 1;
      pows.pows[1] = 0;
      dstCoords[1].AddTerm(pows, 1);

      return;
      //      return Space::Vector2(srcSurfaceCoords.y, srcSurfaceCoords.x);
    }
    break;
    case 1:
    {
      pows.pows[0] = 0;
      pows.pows[1] = 0;
      dstCoords[0].AddTerm(pows, 1);
      pows.pows[0] = 1;
      pows.pows[1] = 0;
      dstCoords[0].AddTerm(pows, -1);
      pows.pows[0] = 0;
      pows.pows[1] = 1;
      dstCoords[0].AddTerm(pows, -1);

      pows.pows[0] = 0;
      pows.pows[1] = 1;
      dstCoords[1].AddTerm(pows, 1);

      return;

      //return Space::Vector2(Scalar(1.0) - srcSurfaceCoords.x - srcSurfaceCoords.y, srcSurfaceCoords.y);
    }
    break;
    case 2:
    {
      pows.pows[0] = 1;
      pows.pows[1] = 0;
      dstCoords[0].AddTerm(pows, 1);

      pows.pows[0] = 0;
      pows.pows[1] = 0;
      dstCoords[1].AddTerm(pows, 1);
      pows.pows[0] = 1;
      pows.pows[1] = 0;
      dstCoords[1].AddTerm(pows, -1);
      pows.pows[0] = 0;
      pows.pows[1] = 1;
      dstCoords[1].AddTerm(pows, -1);

      return;
      //      return Space::Vector2(srcSurfaceCoords.x, Scalar(1.0) - srcSurfaceCoords.x - srcSurfaceCoords.y);
    }
    break;
    }
    assert(0);
  }

  static Scalar GetFaceRefArea(IndexType refFaceNumber)
  {
    Vector p[Space3::NodesPerFace];
    p[0] = GetSurfaceRefCoords(refFaceNumber, Vector2(0.0, 0.0));
    p[1] = GetSurfaceRefCoords(refFaceNumber, Vector2(0.0, 1.0));
    p[2] = GetSurfaceRefCoords(refFaceNumber, Vector2(1.0, 0.0));
    Scalar area = ((p[2] - p[0]) ^ (p[1] - p[0])).Len();
    return fabs(area) * Scalar(0.5);
  }

  static Scalar GetRefCellVolume()
  {
    return Scalar(1.0) / Scalar(6.0);
  }
};

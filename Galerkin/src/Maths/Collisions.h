#pragma once
#include "Spaces.h"

bool CellsCollide(
  const Space2::Vector points0[Space2::NodesPerCell],
  const Space2::Vector points1[Space2::NodesPerCell],
  Space2::Scalar collisionWidth)
{
  Space2::Vector axes[2 * Space2::EdgesPerCell];
  for (int pointNumber = 0; pointNumber < Space2::NodesPerCell; pointNumber++)
  {
    axes[pointNumber] = (points0[(pointNumber + 1) % Space2::NodesPerCell] - points0[pointNumber]).GetPerpendicular().GetNorm();
    axes[pointNumber + Space2::NodesPerCell] = (points1[(pointNumber + 1) % Space2::NodesPerCell] - points1[pointNumber]).GetPerpendicular().GetNorm();
  }

  for (int axisNumber = 0; axisNumber < 2 * Space2::EdgesPerCell; axisNumber++)
  {
    Space2::Scalar mins[2], maxes[2];
    mins[0] = mins[1] = std::numeric_limits<Space2::Scalar>::max() / 2;
    maxes[0] = maxes[1] = -std::numeric_limits<Space2::Scalar>::max() / 2;

    for (int pointNumber = 0; pointNumber < 3; pointNumber++)
    {
      Space2::Scalar proj[2];
      proj[0] = axes[axisNumber] * points0[pointNumber];
      proj[1] = axes[axisNumber] * points1[pointNumber];

      mins[0] = Min(mins[0], proj[0]);
      mins[1] = Min(mins[1], proj[1]);

      maxes[0] = Max(maxes[0], proj[0]);
      maxes[1] = Max(maxes[1], proj[1]);
    }
    if (mins[0] > maxes[1] + collisionWidth || mins[1] > maxes[0] + collisionWidth)
    {
      return false;
    }
  }
  return true;
}

bool CellsCollide(
  const Space3::Vector points0[Space3::NodesPerCell],
  const Space3::Vector points1[Space3::NodesPerCell],
  Space3::Scalar collisionWidth)
{
  Space3::Vector axes[2 * Space3::FacesPerCell + Space3::EdgesPerCell * Space3::EdgesPerCell];
  
  for (int pointNumber = 0; pointNumber < Space3::FacesPerCell; pointNumber++)
  {
    int faceNodes[Space3::NodesPerFace];
    for (int i = 0; i < Space3::NodesPerFace; ++i)
    {
      faceNodes[i] = (pointNumber + i + 1) % Space3::NodesPerCell;
    }
    axes[pointNumber                       ] = ((points0[faceNodes[2]] - points0[faceNodes[0]]) ^ (points0[faceNodes[1]] - points0[faceNodes[0]])).GetNorm();
    axes[pointNumber + Space3::FacesPerCell] = ((points1[faceNodes[2]] - points1[faceNodes[0]]) ^ (points1[faceNodes[1]] - points1[faceNodes[0]])).GetNorm();
  }

  int axesCount = 2 * Space3::FacesPerCell;

  for (int p0 = 0; p0 + 1 < Space3::NodesPerCell; ++p0)
    for (int p1 = p0 + 1; p1 < Space3::NodesPerCell; ++p1)
    {
      Space3::Vector v0 = points0[p1] - points0[p0];
      for (int i0 = 0; i0 + 1 < Space3::NodesPerCell; ++i0)
        for (int i1 = i0 + 1; i1 < Space3::NodesPerCell; ++i1)
        {
          Space3::Vector v1 = points1[i1] - points1[i0];
          axes[axesCount] = (v0 ^ v1).GetNorm();
          axesCount++;
        }
    }

  for (int axisNumber = 0; axisNumber < axesCount; axisNumber++)
  {
    Space3::Scalar mins[2], maxes[2];
    mins[0]  = mins[1]  =  std::numeric_limits<Space3::Scalar>::max() / 2;
    maxes[0] = maxes[1] = -std::numeric_limits<Space3::Scalar>::max() / 2;

    for (int pointNumber = 0; pointNumber < Space3::NodesPerCell; pointNumber++)
    {
      Space3::Scalar proj[2];
      proj[0] = axes[axisNumber] * points0[pointNumber];
      proj[1] = axes[axisNumber] * points1[pointNumber];

      mins[0] = Min(mins[0], proj[0]);
      mins[1] = Min(mins[1], proj[1]);

      maxes[0] = Max(maxes[0], proj[0]);
      maxes[1] = Max(maxes[1], proj[1]);
    }
    if (mins[0] > maxes[1] + collisionWidth || mins[1] > maxes[0] + collisionWidth)
    {
      return false;
    }
  }
  return true;
}

/*
*
*  Two dimensional Triangle-Triangle Overlap Test
*
*/

/* some 2D macros */

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))

#define INTERSECTION_TEST_VERTEXA(P1, Q1, R1, P2, Q2, R2) {\
if (ORIENT_2D(R2,P2,Q1) >= 0.0)\
if (ORIENT_2D(R2,Q2,Q1) <= 0.0)\
if (ORIENT_2D(P1,P2,Q1) > 0.0) {\
if (ORIENT_2D(P1,Q2,Q1) <= 0.0) return 1; \
else return 0;} else {\
if (ORIENT_2D(P1,P2,R1) >= 0.0)\
if (ORIENT_2D(Q1,R1,P2) >= 0.0) return 1; \
else return 0;\
else return 0;}\
else \
if (ORIENT_2D(P1,Q2,Q1) <= 0.0)\
if (ORIENT_2D(R2,Q2,R1) <= 0.0)\
if (ORIENT_2D(Q1,R1,Q2) >= 0.0) return 1; \
else return 0;\
else return 0;\
else return 0;\
else\
if (ORIENT_2D(R2,P2,R1) >= 0.0) \
if (ORIENT_2D(Q1,R1,R2) >= 0.0)\
if (ORIENT_2D(P1,P2,R1) >= 0.0) return 1;\
else return 0;\
else \
if (ORIENT_2D(Q1,R1,Q2) >= 0.0) {\
if (ORIENT_2D(R2,R1,Q2) >= 0.0) return 1; \
else return 0; }\
else return 0; \
else  return 0; \
};

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
  if (ORIENT_2D(R2,P2,Q1) >= 0.0)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0)\
      if (ORIENT_2D(P1,P2,Q1) > 0.0) {\
        if (ORIENT_2D(P1,Q2,Q1) <= 0.0) return 1; \
        else return 0;} else {\
        if (ORIENT_2D(P1,P2,R1) >= 0.0)\
          if (ORIENT_2D(Q1,R1,P2) >= 0.0) return 1; \
          else return 0;\
        else return 0;}\
    else \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0)\
        if (ORIENT_2D(R2,Q2,R1) <= 0.0)\
          if (ORIENT_2D(Q1,R1,Q2) >= 0.0) return 1; \
          else return 0;\
        else return 0;\
      else return 0;\
  else\
    if (ORIENT_2D(R2,P2,R1) >= 0.0) \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0)\
        if (ORIENT_2D(P1,P2,R1) >= 0.0) return 1;\
        else return 0;\
      else \
        if (ORIENT_2D(Q1,R1,Q2) >= 0.0) {\
          if (ORIENT_2D(R2,R1,Q2) >= 0.0) return 1; \
          else return 0; }\
        else return 0; \
    else  return 0; \
 };

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
if (ORIENT_2D(R2,P2,Q1) >= 0.0) {\
if (ORIENT_2D(P1,P2,Q1) >= 0.0) { \
if (ORIENT_2D(P1,Q1,R2) >= 0.0) return 1; \
else return 0;} else { \
if (ORIENT_2D(Q1,R1,P2) >= 0.0){ \
if (ORIENT_2D(R1,P1,P2) >= 0.0) return 1; else return 0;} \
else return 0; } \
} else {\
if (ORIENT_2D(R2,P2,R1) >= 0.0) {\
if (ORIENT_2D(P1,P2,R1) >= 0.0) {\
if (ORIENT_2D(P1,R1,R2) >= 0.0) return 1;  \
else {\
if (ORIENT_2D(Q1,R1,R2) >= 0.0) return 1; else return 0;}}\
else  return 0; }\
else return 0; }}

template <typename T>
int ccw_tri_tri_intersection_2d(const T p1[2], const T q1[2], const T r1[2],
  const T p2[2], const T q2[2], const T r2[2]) {
  if (ORIENT_2D(p2, q2, p1) >= T(0.0)) {
    if (ORIENT_2D(q2, r2, p1) >= T(0.0)) {
      if (ORIENT_2D(r2, p2, p1) >= T(0.0)) return 1;
      else INTERSECTION_TEST_EDGE(p1, q1, r1, p2, q2, r2)
    }
    else {
      if (ORIENT_2D(r2, p2, p1) >= T(0.0))
        INTERSECTION_TEST_EDGE(p1, q1, r1, r2, p2, q2)
      else INTERSECTION_TEST_VERTEX(p1, q1, r1, p2, q2, r2)
    }
  }
  else {
    if (ORIENT_2D(q2, r2, p1) >= T(0.0)) {
      if (ORIENT_2D(r2, p2, p1) >= T(0.0))
        INTERSECTION_TEST_EDGE(p1, q1, r1, q2, r2, p2)
      else  INTERSECTION_TEST_VERTEX(p1, q1, r1, q2, r2, p2)
    }
    else INTERSECTION_TEST_VERTEX(p1, q1, r1, r2, p2, q2)
  }
};

template <typename T>
int tri_tri_overlap_test_2d(const T p1[2], const T q1[2], const T r1[2],
  const T p2[2], const T q2[2], const T r2[2]) {
  if (ORIENT_2D(p1, q1, r1) < T(0.0))
    if (ORIENT_2D(p2, q2, r2) < T(0.0))
      return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, r2, q2);
    else
      return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, q2, r2);
  else
    if (ORIENT_2D(p2, q2, r2) < T(0.0))
      return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, r2, q2);
    else
      return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, q2, r2);
};

bool CellsCollide(
  const Space2::Vector points0[Space2::NodesPerCell],
  const Space2::Vector points1[Space2::NodesPerCell])
{
  return tri_tri_overlap_test_2d<Space2::Scalar>(&(points0[0].x), &(points0[1].x), &(points0[2].x),
    &(points1[0].x), &(points1[1].x), &(points1[2].x)) != 0;
}

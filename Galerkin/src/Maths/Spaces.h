#pragma once

#include "VectorMaths.h"

struct SpaceBase
{
  typedef double             Scalar;
  typedef size_t             IndexType;
};

struct Space1: public SpaceBase
{
  typedef SpaceBase::Scalar             Scalar;
  typedef SpaceBase::IndexType          IndexType;

  const   static IndexType              Dimension = 1;
  typedef Scalar                        Vector;
};

struct Space2: public SpaceBase
{
  typedef SpaceBase::Scalar             Scalar;
  typedef SpaceBase::IndexType          IndexType;

  const   static IndexType              Dimension = 2;
  const   static IndexType              NodesPerCell = 3;
  const   static IndexType              EdgesPerCell = NodesPerCell;
  const   static IndexType              NodesPerEdge = 2;
  const   static IndexType              NodesPerFace = NodesPerEdge; // for uniformity with 3d
  const   static IndexType              FacesPerCell = EdgesPerCell;
  typedef Vector2<IndexType>            IndexVector;
  typedef Vector2<Scalar>               Vector;
  typedef Vector2<Scalar>               Vector2T;
  typedef Vector3<Scalar>               Vector3T;
  typedef Coordinates2<Scalar>          Coords;
  typedef AABB2<Scalar>                 AABB;
  typedef AABB2<IndexType>              IndexAABB;
  typedef AABB2<int>                    IntAABB;
  typedef Tensor2<Scalar>               Tensor;
  typedef Space1                        BorderSpace;

  typedef Vector2 < Vector2<Scalar>>    AsymmetricTensor;

  typedef Scalar                        Angle;
  typedef Scalar                        DiagTensor;
};

struct Space3: public SpaceBase
{
  typedef SpaceBase::Scalar             Scalar;
  typedef SpaceBase::IndexType          IndexType;

  const   static IndexType              Dimension = 3;
  const   static IndexType              NodesPerCell = 4; // dimension + 1
  const   static IndexType              FacesPerCell = NodesPerCell;
  const   static IndexType              EdgesPerCell = 6;
  const   static IndexType              NodesPerFace = 3;
  const   static IndexType              NodesPerEdge = 2;
  typedef Vector3<IndexType>            IndexVector;
  typedef Vector3<Scalar>               Vector;
  typedef Vector2<Scalar>               Vector2T;
  typedef Vector3<Scalar>               Vector3T;
  typedef Coordinates3<Scalar>          Coords;
  typedef AABB3<Scalar>                 AABB;
  typedef AABB3<IndexType>              IndexAABB;
  typedef AABB3<int>                    IntAABB;
  typedef Tensor3<Scalar>               Tensor;
  typedef Vector3<Scalar>               Angle;
  typedef Vector3<Scalar>               DiagTensor;

  typedef Vector3 < Vector3<Scalar> >   AsymmetricTensor;

  typedef Space2                        BorderSpace;
};

template<typename Space>
struct Overload {};

#define SPACE_TYPEDEFS \
  typedef typename Space::Vector        Vector; \
  typedef typename Space::Coords        Coords; \
  typedef typename Space::AABB          AABB; \
  typedef typename Space::Angle         Angle; \
  typedef typename Space::DiagTensor    DiagTensor; \
  typedef typename Space::Tensor        Tensor; \
  typedef typename Space::Scalar        Scalar; \
  typedef typename Space::IndexType     IndexType; \
  typedef typename Space::IndexVector   IndexVector; \
  typedef typename Space::IndexAABB     IndexAABB; \
  typedef typename Space::IntAABB       IntAABB; \
  typedef typename Space::Vector2T      Vector2; \
  typedef typename Space::Vector3T      Vector3; \
  typedef typename Space::BorderSpace   BorderSpace; \
  typedef typename Space::AsymmetricTensor AsymmetricTensor; 

#define SPACE2_TYPEDEFS \
  typedef Space2::Vector        Vector; \
  typedef Space2::Coords        Coords; \
  typedef Space2::AABB          AABB; \
  typedef Space2::Angle         Angle; \
  typedef Space2::DiagTensor    DiagTensor; \
  typedef Space2::Tensor        Tensor; \
  typedef Space2::Scalar        Scalar; \
  typedef Space2::IndexType     IndexType; \
  typedef Space2::IndexVector   IndexVector; \
  typedef Space2::IndexAABB     IndexAABB; \
  typedef Space2::IntAABB       IntAABB; \
  typedef Space2::Vector2T      Vector2; \
  typedef Space2::Vector3T      Vector3; \
  typedef Space2::BorderSpace   BorderSpace; \
  typedef Space2::AsymmetricTensor AsymmetricTensor; 

#define SPACE3_TYPEDEFS \
  typedef Space3::Vector        Vector; \
  typedef Space3::Coords        Coords; \
  typedef Space3::AABB          AABB; \
  typedef Space3::Angle         Angle; \
  typedef Space3::DiagTensor    DiagTensor; \
  typedef Space3::Tensor        Tensor; \
  typedef Space3::Scalar        Scalar; \
  typedef Space3::IndexType     IndexType; \
  typedef Space3::IndexVector   IndexVector; \
  typedef Space3::IndexAABB     IndexAABB; \
  typedef Space3::IntAABB       IntAABB; \
  typedef Space3::Vector2T      Vector2; \
  typedef Space3::Vector3T      Vector3; \
  typedef Space3::BorderSpace   BorderSpace; \
  typedef Space3::AsymmetricTensor AsymmetricTensor; 

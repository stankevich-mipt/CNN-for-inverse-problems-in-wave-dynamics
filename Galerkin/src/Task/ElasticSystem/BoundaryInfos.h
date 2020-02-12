#pragma once
#include "../../Maths/Spaces.h"

// boundaries
struct BoundaryConditions
{
  enum Types
  {
    Absorb = 0, Free = 1, Fixed, Symmetry, AntiSymmetry
  };
};

struct BoundaryDescription
{
  BoundaryDescription(): infoIndex(-1),
    reflectionCoeff(1.0), type(BoundaryConditions::Absorb)
  {}

  int infoIndex;
  double reflectionCoeff;
  BoundaryConditions::Types type;
};

template <typename Space>
struct BoundaryInfoFunctor
{
  SPACE_TYPEDEFS

  virtual ~BoundaryInfoFunctor() = default;
  virtual void operator()(const Vector& globalPoint,
                          const Vector& externalNormal,
                          const Scalar currTime,
                          Scalar* values) = 0;

  virtual void SetCurrentVelocity(const Vector& velocity) {}
};

template <typename Space>
struct FreeBoundaryInfo
{
  SPACE_TYPEDEFS
  FreeBoundaryInfo(): 
    boundaryFunctor(nullptr), 
    dynamicContactInteractionType(IndexType(-1))
  {
  }
  BoundaryInfoFunctor<Space>* boundaryFunctor;
  IndexType dynamicContactInteractionType;
};

template <typename Space>
struct FixedBoundaryInfo
{
  SPACE_TYPEDEFS

  FixedBoundaryInfo(): boundaryFunctor(nullptr)
  {
  }
  BoundaryInfoFunctor<Space>* boundaryFunctor;
};

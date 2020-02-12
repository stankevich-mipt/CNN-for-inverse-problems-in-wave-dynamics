#pragma once

#include "../../Maths/Spaces.h"

// contacts
struct ContactConditions
{
  enum Types
  {
    Glue = 0, Glide = 1, Friction = 2
  };
};

struct ContactDescription
{
  ContactDescription(): 
    type(ContactConditions::Glue), infoIndex(-1)
  {}

  ContactConditions::Types type;
  int infoIndex;
};

template <typename Space>
struct GlueContactInfo
{
  SPACE_TYPEDEFS

  GlueContactInfo(): dynamicBoundaryInteractionType(IndexType(-1))
  {
  }

  GlueContactInfo(Scalar maxShearStress, Scalar maxLongitudinalStress, 
    IndexType dynamicBoundaryInteractionType):
    maxShearStress(maxShearStress),
    maxLongitudinalStress(maxLongitudinalStress),
    dynamicBoundaryInteractionType(dynamicBoundaryInteractionType)
  {
  }
  Scalar maxShearStress; 
  Scalar maxLongitudinalStress;
  IndexType dynamicBoundaryInteractionType;
};

template <typename Space>
struct FrictionContactInfo
{
  SPACE_TYPEDEFS
  FrictionContactInfo()
  {
    frictionCoeff = 0;
    dynamicBoundaryInteractionType = -1;
  }
  FrictionContactInfo(Scalar frictionCoeff, IndexType dynamicBoundaryInteractionType):
    frictionCoeff(frictionCoeff),
    dynamicBoundaryInteractionType(dynamicBoundaryInteractionType)
  {
  }
  Scalar frictionCoeff; 
  IndexType dynamicBoundaryInteractionType;
};


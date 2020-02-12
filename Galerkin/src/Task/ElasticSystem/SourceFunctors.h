#pragma once
#include <limits>

#include "../../Maths/Spaces.h"

template <typename Space>
struct SourceFunctor
{
  SPACE_TYPEDEFS

  SourceFunctor(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult):
    tensionDimensionlessMult(tensionDimensionlessMult),
    velocityDimensionlessMult(velocityDimensionlessMult)
  {
  }

  virtual ~SourceFunctor() = default;
  virtual void operator ()(const Vector& point, Scalar time, Scalar* values) const = 0;

  Scalar tensionDimensionlessMult;
  Scalar velocityDimensionlessMult;
};

template <typename ElasticSystem>
struct ConstantAccelerationFunctor: public SourceFunctor<typename ElasticSystem::Space>
{
  typedef typename ElasticSystem::Space Space;
  SPACE_TYPEDEFS

  const static IndexType dimsCount = ElasticSystem::dimsCount;
  
  using SourceFunctor<Space>::tensionDimensionlessMult;
  using SourceFunctor<Space>::velocityDimensionlessMult;

  ConstantAccelerationFunctor(const Vector& acceleration, Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult): 
    SourceFunctor<typename ElasticSystem::Space>(tensionDimensionlessMult, velocityDimensionlessMult), 
    acceleration(acceleration)
  {
  }

  virtual ~ConstantAccelerationFunctor() = default;
  virtual void operator()(const Vector&, Scalar, Scalar* values) const override
  {
    std::fill_n(values, dimsCount, Scalar(0.0));
    SetAcceleration(Overload<Space>(), values);
  }

protected: 
  void SetAcceleration(Overload<Space2>, typename Space2::Scalar* values) const
  {
    values[3] = acceleration.x / this->velocityDimensionlessMult;
    values[4] = acceleration.y / this->velocityDimensionlessMult;
  }

  void SetAcceleration(Overload<Space3>, typename Space3::Scalar* values) const
  {
    values[6] = acceleration.x / this->velocityDimensionlessMult;
    values[7] = acceleration.y / this->velocityDimensionlessMult;
    values[8] = acceleration.z / this->velocityDimensionlessMult;
  }

  Vector acceleration;
};
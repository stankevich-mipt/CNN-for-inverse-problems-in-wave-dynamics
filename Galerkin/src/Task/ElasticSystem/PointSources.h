#pragma once
#include <limits>

template <typename Space>
class PointSource
{
public:
  SPACE_TYPEDEFS

  PointSource(Vector point, Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult):
    point(point),
    tensionDimensionlessMult(tensionDimensionlessMult),
    velocityDimensionlessMult(velocityDimensionlessMult)
  {}

  virtual ~PointSource() = default;
  virtual Vector GetPoint() const
  {
    return point;
  }

  virtual void operator()(Scalar time, Scalar* values) const = 0;

protected:
  Vector point;
  Scalar tensionDimensionlessMult;
  Scalar velocityDimensionlessMult;
};

template <typename ElasticSystem>
class ForcePointSource: public PointSource<typename ElasticSystem::Space>
{
public:
  typedef typename ElasticSystem::Space Space;
  SPACE_TYPEDEFS
  const static IndexType dimsCount = ElasticSystem::dimsCount;
  using PointSource<typename ElasticSystem::Space>::point;

  ForcePointSource(const Vector& point,
    Scalar peakFrequency, Vector acceleration, Scalar latency,
    Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult): 
    PointSource<Space>(point, tensionDimensionlessMult, velocityDimensionlessMult), 
    peakFrequency(peakFrequency), acceleration(acceleration), latency(latency)
  {
    this->latency = std::max(latency, sqrt(1.5) / (pi * peakFrequency) * 3);
  }

  void operator()(Scalar time, Scalar* values) const override
  {
    Scalar arg = Sqr(pi * peakFrequency * (time - latency));
    Scalar mult = (1 - 2 * arg) * exp(-arg);
    std::fill_n(values, dimsCount, Scalar(0.0));
    
    SetVelocity(Overload<Space>(), mult, values);
  }

private:
  void SetVelocity(Overload<Space2>, typename Space2::Scalar mult, typename Space2::Scalar* values) const
  {
    values[3] = mult * acceleration.x;
    values[4] = mult * acceleration.y;
  }

  void SetVelocity(Overload<Space3>, typename Space3::Scalar mult, typename Space3::Scalar* values) const
  {
    values[6] = mult * acceleration.x;
    values[7] = mult * acceleration.y;
    values[8] = mult * acceleration.z;
  }

  Scalar peakFrequency;
  Vector acceleration;
  Scalar latency; 
};

// spherical point explosion
template <typename ElasticSystem>
struct MonopoleSource: public PointSource<typename ElasticSystem::Space>
{
  typedef typename ElasticSystem::Space Space;
  SPACE_TYPEDEFS

  const static IndexType dimsCount = ElasticSystem::dimsCount;

  using PointSource<Space>::tensionDimensionlessMult;
  using PointSource<Space>::velocityDimensionlessMult;
  using PointSource<Space>::point;

  MonopoleSource(const Vector& point, Scalar pressure, Scalar peakFrequency, Scalar latency, 
    Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult):
    PointSource<Space>(point, tensionDimensionlessMult, velocityDimensionlessMult),
    pressure(pressure),
    peakFrequency(peakFrequency),
    latency(latency)
  {
    this->latency = std::max(latency, sqrt(1.5) / (pi * peakFrequency) * 3);
  }

  virtual ~MonopoleSource() = default;

  virtual void operator()(Scalar time, Scalar* values) const override
  {
    std::fill_n(values, dimsCount, Scalar(0.0));

    Scalar arg = Sqr(pi * peakFrequency * (time - latency));
    Scalar mult = (1 - 2 * arg) * exp(-arg);
    Scalar p = pressure * mult;

    for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
    {
      values[dimIndex] = -p / Space::Dimension; // sigma_ii
    }
  }

protected:
  Scalar pressure;
  Scalar peakFrequency;
  Scalar latency;
};



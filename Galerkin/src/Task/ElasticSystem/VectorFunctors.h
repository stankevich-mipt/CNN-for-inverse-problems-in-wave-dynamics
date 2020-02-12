#pragma once

#include <limits>
#include <vector>
#include "../../Maths/Spaces.h"

template <typename Space>
struct VectorFunctor
{
  SPACE_TYPEDEFS

  virtual ~VectorFunctor() = default;
  virtual Vector operator ()(const Vector& point, const Vector& norm, Scalar time) const = 0;

  virtual void SetCurrentVelocity(const Vector& velocity) {}
};

template <typename Space>
struct ConstVectorFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS
  
  ConstVectorFunctor(Vector constantValue)
  {
    this->constantValue = constantValue;
  }

  Vector operator()(const Vector&, const Vector&, Scalar) const override
  {
    return constantValue;
  }

protected:
  Vector constantValue;
};

template <typename Space>
struct BoxVectorFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  BoxVectorFunctor(AABB aabb, Vector constantValue, Vector aabbVelocity)
  {
    this->constantValue = constantValue;
    this->aabb          = aabb;
    this->aabbVelocity  = aabbVelocity;
  }

  Vector operator()(const Vector& point, const Vector&, Scalar time) const override
  {
    AABB movedAABB;
    movedAABB.Expand(aabb.boxPoint1 + aabbVelocity * time);
    movedAABB.Expand(aabb.boxPoint2 + aabbVelocity * time);
    if (movedAABB.Includes(point))
      return constantValue;
    else
      return Space::Vector::zeroVector();
  }

private:
  Vector constantValue;
  AABB   aabb;
  Vector aabbVelocity;
};

template <typename Space>
struct RadialVectorFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  RadialVectorFunctor(Vector center, Scalar intensity)
  {
    this->center    = center;
    this->intensity = intensity;
  }

  Vector operator()(const Vector& point, const Vector&, Scalar) const override
  {
    return (point - center).GetNorm() * intensity;
  }

private:
  Vector center;
  Scalar intensity;
};

template <typename Space>
struct WaveVectorFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  WaveVectorFunctor(const Vector& waveVector, Scalar waveLength, Scalar initialPhase, Scalar speed, 
                    Scalar maxTime = std::numeric_limits<Scalar>::infinity(), bool shear = false):
    invWaveLength(Scalar(1.0) / waveLength),
    initialPhase(initialPhase),
    speed(speed),
    maxTime(maxTime),
    shear(shear)
  {
    normal = waveVector.GetNorm();
    magnitude = waveVector.Len();
    tangent = waveVector.GetNorm().GetPerpendicular();
    k = 2 * Scalar(pi) * invWaveLength;
    omega = 2 * Scalar(pi) * speed * invWaveLength;
  }

  Vector operator()(const Vector& point, const Vector& faceNormal, Scalar time) const override
  {
    Vector value;

    if (time < maxTime || maxTime < 0)
    {
      Scalar mult = sin(point * normal * k - omega * time + initialPhase);
      Scalar modulatedMagnitude = magnitude /** (1 - fabs(2 * time - maxTime) / maxTime)*/; //saw
      if (shear)
      {
        value = tangent * modulatedMagnitude * mult;
      } else
      {     
      // normal -> faceNormal
        value = normal * modulatedMagnitude * mult;
      }
    } else
    {
      value = Space::Vector::zeroVector();
    }
    return value;
  }

private:
  Vector normal;
  Vector tangent;
  Scalar invWaveLength;
  Scalar initialPhase;
  Scalar speed;
  Scalar magnitude;
  Scalar maxTime;
  bool shear;
  Scalar k;
  Scalar omega;
};

template <typename Space>
struct RotatingVectorFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  RotatingVectorFunctor(const Vector& pos, Angle angularVelocity, Vector linearVelocity, Scalar linearTime):
    pos(pos),
    angularVelocity(angularVelocity),
    linearVelocity(linearVelocity),
    linearTime(linearTime)
  {
  }

  Vector operator()(const Vector& point, const Vector& faceNormal, Scalar time) const override
  {
    Vector value;

    Vector currLinearVelocity = Vector::zero();
    Vector displacement       = Vector::zero();

    if (linearTime < 0 || time < linearTime)
    {
      currLinearVelocity = linearVelocity;
      displacement = linearVelocity * time;
    } else
    {
      displacement = linearVelocity * std::max(linearTime, Scalar(0));
    }

    // TODO: fix for 3d case
    // Vector tangent = (point - (pos + displacement)).GetPerpendicular();
    // Vector velocity = currLinearVelocity + tangent * angularVelocity;

    Vector velocity = Vector::zero();
    return velocity;
  }

private:
  Vector pos;
  Angle angularVelocity;
  Vector linearVelocity;
  Scalar linearTime;
};

template <typename Space>
struct RiekerVectorFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  // See alse: http://subsurfwiki.org/wiki/Ricker_wavelet
  RiekerVectorFunctor(const Vector& waveVector, Scalar peakFrequency):
    waveVector(waveVector), peakFrequency(peakFrequency)
  {
    delay = sqrt(1.5) / (pi * peakFrequency) * 3;
  }

  Vector operator()(const Vector&, const Vector&, Scalar time) const override
  {
    Scalar arg = Sqr(pi * peakFrequency * (time - delay));
    Scalar mult = (1 - 2 * arg) * exp(-arg);
    return waveVector * mult;
  }

  Vector waveVector;
  Scalar  peakFrequency;
  Scalar  delay;
};


template <typename Space>
struct TabulatedVectorFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  struct Sample
  {
    Scalar time;
    Scalar value;
    bool operator<(const Sample& other) const
    {
      return time < other.time;
    }
  };

  TabulatedVectorFunctor(const Vector& direction, const std::string& fileName):
    direction(direction), fileName(fileName)
  {
    std::fstream input(fileName.c_str());
    assert(input.fail() == false);

    this->direction.Normalize();

    IndexType samplesCount;
    input >> samplesCount;
    samples.reserve(samplesCount);
    for (IndexType sampleIndex = 0; sampleIndex < samplesCount; ++sampleIndex)
    {
      Sample sample;
      input >> sample.time >> sample.value;
      samples.push_back(sample);
    }
    std::sort(samples.begin(), samples.end());
  }

  Vector operator()(const Vector&, const Vector&, Scalar time) const override
  {
    IndexType sampleIndex = 0;
    IndexType leftBorder = 0;
    IndexType rightBorder = samples.size();
    while (leftBorder + 1 < rightBorder)
    {
      sampleIndex = (leftBorder + rightBorder) / 2;
      samples[sampleIndex].time > time ? rightBorder : leftBorder = sampleIndex;
    }
    if (sampleIndex == samples.size()) return Vector::zero();
    Scalar mult = (samples[sampleIndex + 1].value - samples[sampleIndex].value) /
      (samples[sampleIndex + 1].time - samples[sampleIndex].time) * (time - samples[sampleIndex].time) + samples[sampleIndex].value;
    return direction * mult;
  }

  std::vector<Sample> samples;
  Vector direction;
  std::string fileName;
};

template <typename Space>
struct HydraulicPressureFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  HydraulicPressureFunctor(Scalar fluidRho, Vector g, Vector fluidSurfacePoint):
    fluidRho(fluidRho), g(g), fluidSurfacePoint(fluidSurfacePoint)
  {}

  Vector operator()(const Vector& point, const Vector& normal, Scalar) const override
  {
    Scalar pressure = fluidRho * g * (point - fluidSurfacePoint);
    if (pressure > 0)
    {
      return -normal * pressure;
    } else
    {
      return Vector::zero();
    }
  }

  Scalar  fluidRho;
  Vector g;
  Vector fluidSurfacePoint;
};

template <typename Space>
struct HydrodynamicResistanceFunctor: public VectorFunctor<Space>
{
  SPACE_TYPEDEFS

  HydrodynamicResistanceFunctor(AABB aabb, Scalar mediumRho, Vector flowVelocity, Scalar alpha, Scalar beta):
    aabb(aabb), mediumRho(mediumRho), flowVelocity(flowVelocity), alpha(alpha), beta(beta)
  {}
  
  virtual void SetCurrentVelocity(const Vector& velocity) override
  {
    this->v = velocity;
  }

  Vector operator()(const Vector& point, const Vector& n, Scalar /* time */) const override
  { 
    Scalar vn = (v - flowVelocity) * n;
    if (aabb.Includes(point) && vn >= 0)
    {
      Vector tau = n.GetPerpendicular();
      Scalar vtau = (v - flowVelocity) * tau;
      Vector res = -Scalar(0.5) * mediumRho * (alpha * Sqr(vn) * n + beta * Sqr(vtau) * Sgn(vtau) * tau);
      return res;
    } else
    {
      return Vector::zero();
    }
  }

  AABB    aabb;
  Scalar  mediumRho;
  Vector  flowVelocity;
  Scalar  alpha;
  Scalar  beta;
  Vector  v;
};

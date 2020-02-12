#pragma once

#include <limits>
#include <stdexcept>
#include "../../Maths/Spaces.h"

template<typename ElasticSpace>
class IniStateMaker
{
public:
  typedef typename ElasticSpace::SpaceType Space;
  SPACE_TYPEDEFS

  virtual ~IniStateMaker() = default;
  virtual void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) = 0;

  // point, lambda, mju, invRho have to be dimensionless
  virtual typename ElasticSpace::Elastic GetValue(const Vector& point,
    const Scalar lambda, const Scalar mju, const Scalar invRho) = 0;
};


template<typename ElasticSpace>
class RickerIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  SPACE_TYPEDEFS

  RickerIniStateMaker(Vector pos, Vector velocity, Scalar wavelength)
  {
    this->pos = pos;
    this->velocity = velocity;
    this->wavelength = wavelength;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->pos /= velocityDimensionlessMult;
    this->velocity /= velocityDimensionlessMult;
    this->wavelength /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Scalar pSpeed = sqrt((lambda + Scalar(2.0) * mju) * invRho);
    Scalar sSpeed = sqrt(mju * invRho);

    Scalar rhoc1mc3 = (pSpeed - sSpeed) / invRho;
    Scalar rhoc2 = sSpeed / invRho;
    Scalar rhoc3 = sSpeed / invRho;

    Scalar er = wavelength / pi;

    Scalar absv = velocity.Len();
    Vector n = velocity.GetNorm();

    Scalar r = (point - pos) * n * (Scalar(1.0) / er);
    Scalar ex = (1 - 2 * r * r) * exp(- r * r);

    Elastic res;
    res.SetVelocity(ex * velocity);
    res.SetTension(GetTension(Overload<Space>(), rhoc1mc3, rhoc2, rhoc3, absv, n, ex));
    return res;
  }

  typename Space2::Tensor GetTension(Overload<Space2>,
    Scalar rhoc1mc3, Scalar rhoc2, Scalar rhoc3, Scalar absv, typename Space2::Vector n, Scalar ex)
  {
    Scalar xx = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.x * n.x;
    Scalar yy = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.y * n.y;
    Scalar xy = ex * absv * rhoc1mc3 * n.x * n.y;

    return typename Space2::Tensor(xx, xy, yy);
  }

  typename Space3::Tensor GetTension(Overload<Space3>,
    Scalar rhoc1mc3, Scalar rhoc2, Scalar rhoc3, Scalar absv, typename Space3::Vector n, Scalar ex)
  {
    Scalar xx = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.x * n.x;
    Scalar yy = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.y * n.y;
    Scalar zz = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.z * n.z;
    Scalar xy = - ex * absv * rhoc1mc3 * n.x * n.y;
    Scalar xz = - ex * absv * rhoc1mc3 * n.x * n.z;
    Scalar yz = - ex * absv * rhoc1mc3 * n.y * n.z;

    return typename Space3::Tensor(xx, xy, xz, yy, yz, zz);
  }

private:
  Vector pos;
  Vector velocity;
  Scalar wavelength;
};


// plane wave in acoustic medium
template <typename ElasticSpace>
class AcousticPlaneWaveIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::Vector                    Vector;
  typedef typename ElasticSpace::Tensor                    Tensor;
  typedef typename ElasticSpace::Scalar                    Scalar;
  typedef typename ElasticSpace::IndexType                 IndexType;
  typedef typename ElasticSpace::MediumParametersType      MediumParameters;

  AcousticPlaneWaveIniStateMaker(Vector pos, Vector velocity, Scalar waveLength)
  {
    this->pos = pos;
    this->velocity = velocity;
    this->waveLength = waveLength;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->pos /= velocityDimensionlessMult;
    this->velocity /= velocityDimensionlessMult;
    this->waveLength /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar , const Scalar invRho) override
  {
    MediumParameters mediumParameters(lambda, Scalar(0.0), invRho);
    Elastic          elastic;

    // point to wave front line distanse
    Scalar d = (point - pos) * velocity.GetNorm();

    if (fabs(d) <= waveLength / Scalar(2.0)) 
    {
      Scalar mag = velocity.Len() * GetGauss((waveLength / Scalar(2.0) + d) / waveLength, Scalar(0.5), Scalar(0.2));
      elastic.SetVelocity(velocity.GetNorm() * mag);
      Scalar pressure = mag * mediumParameters.GetZp();
      /*
        pressure > 0 - extension
        pressure < 0 - compression
        In gases and fluids - only compression
      */
      elastic.SetTension(Tensor(-pressure));
    } else
    {
      elastic.SetZeroValues();
    }
    return elastic;
  }

  Scalar GetGauss(Scalar x, Scalar mju, Scalar sigma) const
  {
    return exp(-Scalar(0.5) * Sqr((x - mju) / sigma) / (sigma * sqrt(2 * pi)));
  }

private:
  Vector pos;
  Vector velocity;
  Scalar  waveLength;
};

template<typename ElasticSpace>
class BerlageIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  SPACE_TYPEDEFS

  BerlageIniStateMaker(Vector pos, Vector velocity, Scalar wavelength, Scalar smoothFactor, bool shear)
  {
    this->pos = pos;
    this->velocity = velocity;
    this->wavelength = wavelength;
    this->smoothFactor = smoothFactor;
    this->shear = shear;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->pos /= velocityDimensionlessMult;
    this->velocity /= velocityDimensionlessMult;
    this->wavelength /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Scalar pSpeed = sqrt((lambda + Scalar(2.0) * mju) * invRho);
    Scalar sSpeed = sqrt(mju * invRho);

    Scalar rhoc1mc3 = (pSpeed - sSpeed) / invRho;
    Scalar rhoc2 = sSpeed / invRho;
    Scalar rhoc3 = sSpeed / invRho;
    Scalar er = wavelength / pi;

    Scalar absv = velocity.Len();
    Vector n = velocity.GetNorm();
    Scalar r2 = Sqr((point - pos) * n);
    Scalar ex = 0;

    Scalar a = r2 / Sqr(er);
    ex = (a <= Sqr(1 + smoothFactor)) ? exp(-a) : 0;

    return GetElastic(Overload<Space>(), ex, absv, rhoc1mc3, rhoc2, rhoc3, n);
  }

private:
  Elastic GetElastic(Overload<Space2>, Scalar ex, Scalar absv, 
    Scalar rhoc1mc3, Scalar rhoc2, Scalar rhoc3, const typename Space2::Vector& n)
  {
    Elastic res;
    if (!shear)
    {
      Scalar xx = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.x * n.x;
      Scalar yy = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.y * n.y;
      Scalar xy = -ex * absv * rhoc1mc3 * n.x * n.y;
      res.SetVelocity(ex * velocity);
      res.SetTension(xx, yy, xy);
    } else
    {
      Vector n1 = Vector(n.y, -n.x);//n.GetPerpendicular();
      n1 *= ex * absv / (n1).Len();

      Scalar xx = -rhoc2 * n.x * n1.x * Scalar(2.0);
      Scalar yy = -rhoc2 * n.y * n1.y * Scalar(2.0);
      Scalar xy = -rhoc2 * (n.x * n1.y + n.y * n1.x);
      res.SetVelocity(n1);
      res.SetTension(xx, yy, xy);
    }
    return res;
  }

  Elastic GetElastic(Overload<Space3>, Scalar ex, Scalar absv, 
    Scalar rhoc1mc3, Scalar rhoc2, Scalar rhoc3, const typename Space3::Vector& n)
  {
    Elastic res;
    if (!shear)
    {
      Scalar xx = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.x * n.x;
      Scalar yy = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.y * n.y;
      Scalar zz = -ex * absv * rhoc3 - ex * absv * rhoc1mc3 * n.z * n.z;
      Scalar xy = -ex * absv * rhoc1mc3 * n.x * n.y;
      Scalar xz = -ex * absv * rhoc1mc3 * n.x * n.z;
      Scalar yz = -ex * absv * rhoc1mc3 * n.y * n.z;
      res.SetVelocity(ex * velocity);
      res.SetTension(xx, yy, zz, xy, yz, xz);
    } else
    {
      Vector3 n1 = Vector3(0, -1, 0);//n.GetPerpendicular();
      n1 *= ex * absv / (n1).Len();

      Scalar xx = -rhoc2 * n.x * n1.x * Scalar(2.0);
      Scalar yy = -rhoc2 * n.y * n1.y * Scalar(2.0);
      Scalar zz = -rhoc2 * n.z * n1.z * Scalar(2.0);
      Scalar xy = -rhoc2 * (n.x * n1.y + n.y * n1.x);
      Scalar xz = -rhoc2 * (n.x * n1.z + n.z * n1.x);
      Scalar yz = -rhoc2 * (n.y * n1.z + n.z * n1.y);
      res.SetVelocity(n1);
      res.SetTension(xx, yy, zz, xy, yz, xz);
    }
    return res;
  }

  Vector pos;
  Vector velocity;
  Scalar wavelength;
  Scalar smoothFactor;
  bool shear;
};


template<typename ElasticSpace>
class RadialIniStateMaker : public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  SPACE_TYPEDEFS

  RadialIniStateMaker(Vector pos, Scalar radialVelocity, Scalar wavelength, Scalar smoothFactor)
  {
    this->pos            = pos;
    this->radialVelocity = radialVelocity;
    this->wavelength     = wavelength;
    this->smoothFactor   = smoothFactor;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->pos /= velocityDimensionlessMult;
    this->radialVelocity /= velocityDimensionlessMult;
    this->wavelength /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Scalar pSpeed = sqrt((lambda + Scalar(2.0) * mju) * invRho);
    Scalar sSpeed = sqrt(mju * invRho);

    Scalar rhoc1mc3 = (pSpeed - sSpeed) / invRho;
    Scalar rhoc2 = sSpeed / invRho;
    Scalar rhoc3 = sSpeed / invRho;
    Scalar er = wavelength / pi;

    Vector n;
    if ((point - pos).SquareLen() > 0)
      n = (point - pos).GetNorm();
    else
      n = Vector::xAxis();
    Scalar r2 = (point - pos).SquareLen();
    Scalar ex = 0;

    Scalar a = r2 / Sqr(er);
    ex = (a <= Sqr(1 + smoothFactor)) ? exp(-a) : 0;

    Elastic res;
    res.SetVelocity(ex * radialVelocity * n);
    res.SetTension(GetTension(Overload<Space>(), ex, rhoc1mc3, rhoc2, rhoc3, n));

    return res;
  }
private:
  Tensor GetTension(Overload<Space2>, Scalar ex, Scalar rhoc1mc3, Scalar rhoc2, Scalar rhoc3, const Vector& n)
  {
    Scalar xx = -ex * radialVelocity * rhoc3 - ex * radialVelocity * rhoc1mc3 * n.x * n.x;
    Scalar yy = -ex * radialVelocity * rhoc3 - ex * radialVelocity * rhoc1mc3 * n.y * n.y;
    Scalar xy = -ex * radialVelocity * rhoc1mc3 * n.x * n.y;
    return Tensor(xx, yy, xy);
  }

  typename Space3::Tensor GetTension(Overload<Space3>, Scalar ex, Scalar rhoc1mc3, Scalar rhoc2, Scalar rhoc3, const Space3::Vector& n)
  {
    Scalar xx = -ex * radialVelocity * rhoc3 - ex * radialVelocity * rhoc1mc3 * n.x * n.x;
    Scalar yy = -ex * radialVelocity * rhoc3 - ex * radialVelocity * rhoc1mc3 * n.y * n.y;
    Scalar zz = -ex * radialVelocity * rhoc3 - ex * radialVelocity * rhoc1mc3 * n.z * n.z;
    Scalar xy = -ex * radialVelocity * rhoc1mc3 * n.x * n.y;
    Scalar xz = -ex * radialVelocity * rhoc1mc3 * n.x * n.z;
    Scalar yz = -ex * radialVelocity * rhoc1mc3 * n.y * n.z;
    return typename Space3::Tensor(xx, xy, xz, yy, yz, zz);
  }

  Vector pos;
  Scalar radialVelocity;
  Scalar wavelength;
  Scalar smoothFactor;
};

template<typename ElasticSpace>
class BallIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  SPACE_TYPEDEFS

  BallIniStateMaker(Vector pos, Scalar radialVelocity, Scalar pressure, Scalar wavelength, Scalar smoothFactor)
  {
    this->pos = pos;
    this->pressure = pressure;
    this->radialVelocity = radialVelocity;
    this->wavelength = wavelength;
    this->smoothFactor = smoothFactor;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->pos /= velocityDimensionlessMult;
    this->pressure /= tensionDimensionlessMult;
    this->radialVelocity /= velocityDimensionlessMult;
    this->wavelength /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Scalar pSpeed = sqrt((lambda + Scalar(2.0) * mju) * invRho);
    Scalar sSpeed = sqrt(mju * invRho);

    Scalar rhoc1mc3 = (pSpeed - sSpeed) / invRho;
    Scalar rhoc2 = sSpeed / invRho;
    Scalar rhoc3 = sSpeed / invRho;
    Scalar er = wavelength / pi;

    Vector n;
    if ((point - pos).SquareLen() > 0)
      n = (point - pos).GetNorm();
    else
      n = Vector::xAxis();
    Scalar r2 = (point - pos).SquareLen();
    Scalar ex = 0;
    Scalar a = r2 / Sqr(er);
    ex = (a <= Sqr(1 + smoothFactor)) ? exp(-a) : 0;
    Elastic res;
    res.SetVelocity(ex * n * radialVelocity);
    res.SetTension(Tensor(-ex * pressure));

    return res;
  }

private:
  Vector pos;
  Scalar radialVelocity;
  Scalar pressure;
  Scalar wavelength;
  Scalar smoothFactor;
};

template<typename ElasticSpace>
class BagelIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  typedef typename ElasticSpace::MediumParametersType      MediumParameters;
  SPACE_TYPEDEFS

  BagelIniStateMaker(Scalar r1, Scalar r2, Vector center, Scalar velocityMagnitude)
  {
    if (r2 <= r1) 
    {
      throw std::logic_error("BagelIniStateMaker: r2 <= r1");
    }
    this->r1        = r1;
    this->r2        = r2;
    this->center    = center;
    this->velocityMagnitude = velocityMagnitude;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->r1                /= velocityDimensionlessMult;
    this->r2                /= velocityDimensionlessMult;
    this->center            /= velocityDimensionlessMult;
    this->velocityMagnitude /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    MediumParameters mediumParameters(lambda, mju, invRho);
    Elastic elastic;
    Vector d = point - center;
    Scalar l = d.Len();
    if (l >= r1 && l <= r2)
    {
      Scalar arg = (l - r1) / (r2 - r1);
      Scalar value = velocityMagnitude * GetGauss(arg, Scalar(0.5), Scalar(0.2));
      elastic.SetVelocity(value * d.GetNorm());
      Scalar c1 = mediumParameters.GetPSpeed();
      Scalar c2 = mediumParameters.GetSSpeed();
      Scalar c3 = c1 - 2.0 * c2 * c2 / c1;
      Tensor t(-value * c3 / invRho);
      t += TensorProduct(d.GetNorm()) * (-value / invRho * (c1 - c3));
      elastic.SetTension(t);
    } else
    {
      elastic.SetZeroValues();
    }
    return elastic;
  }

private:
  Scalar r1;
  Scalar r2;
  Scalar velocityMagnitude;
  Vector center;

  Scalar GetGauss(Scalar x, Scalar mju, Scalar sigma) const
  {
    return exp(-Scalar(0.5) * Sqr((x - mju) / sigma) / (sigma * sqrt(2 * pi)));
  }

  Scalar GetRicker(Scalar x, Scalar mju, Scalar sigma) const
  {
    // TODO
    return 1;
  }
};

template<typename ElasticSpace>
class DirectionalExplosionIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  typedef typename ElasticSpace::MediumParametersType      MediumParameters;
  SPACE_TYPEDEFS

  DirectionalExplosionIniStateMaker(Vector pos, Vector velocity, Scalar wavelength, Scalar smoothFactor)
  {
    this->pos = pos;
    this->velocity = velocity;
    this->wavelength = wavelength;
    this->smoothFactor = smoothFactor;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->pos        /= velocityDimensionlessMult;
    this->velocity   /= velocityDimensionlessMult;
    this->wavelength /= velocityDimensionlessMult; 
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Scalar pSpeed = sqrt((lambda + Scalar(2.0) * mju) * invRho);
    Scalar sSpeed = sqrt(mju * invRho);

    Scalar rhoc1mc3 = (pSpeed - sSpeed) / invRho;
    Scalar rhoc2 = sSpeed / invRho;
    Scalar rhoc3 = sSpeed / invRho;
    Scalar er = wavelength / pi;

    Vector n;
    if ((point - pos).SquareLen() > 0)
      n = (point - pos).GetNorm();
    else
      n = Vector::xAxis();
    Scalar r2 = (point - pos).SquareLen();
    Scalar ex = 0;
    Scalar a = r2 / Sqr(er);
    ex = (a <= Sqr(1 + smoothFactor)) ? exp(-a) : 0;
    Elastic res;

    res.SetVelocity(ex * velocity);
    res.SetTension(Tensor(0));

    return res;
  }
private:
  Vector pos;
  Vector velocity;
  Scalar wavelength;
  Scalar smoothFactor;
};

template<typename ElasticSpace>
class BoxIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::SpaceType                 Space;
  SPACE_TYPEDEFS

  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::MediumParametersType      MediumParameters;

  BoxIniStateMaker(Vector velocity, const AABB& fullBoundary)
  {
    this->velocity = velocity;
    this->fullBoundary = fullBoundary;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    this->velocity                        /= velocityDimensionlessMult;
    this->fullBoundary.boxPoint1          /= velocityDimensionlessMult;
    this->fullBoundary.boxPoint2          /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Elastic res;
    res.SetZeroValues();

    if (fullBoundary.Includes(point))
    {
      res.SetVelocity(velocity);
    }
    return res;
  }
private:
  Vector velocity;
  AABB fullBoundary;
};

template<typename ElasticSpace>
class MonochromaticPlaneWaveIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  typedef typename ElasticSpace::MediumParametersType      MediumParameters;
  SPACE_TYPEDEFS

  MonochromaticPlaneWaveIniStateMaker(const Vector& normal, Scalar waveLength, Scalar initialPhase, bool shear = false):
    normal(normal.GetNorm()),
    invWaveLength(Scalar(1.0) / waveLength),
    initialPhase(initialPhase),
    shear(shear)
  {
    tangent = normal.GetPerpendicular();
    k = 2 * pi * invWaveLength;
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    invWaveLength *= velocityDimensionlessMult;
    k             *= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Elastic res; 
    /*
    if (shear)
    {
      Scalar sSpeed = sqrt(mju * invRho);
      Scalar omega = 2 * pi * sSpeed * invWaveLength;

      Scalar mult = cos(point * normal * k - omega * 0 // time 
                        + initialPhase);

      res.SetVelocity(tangent * sSpeed * mult);

      Scalar sigmaxx = -mult * mju * (2 * normal.x * tangent.x);
      Scalar sigmaxy = -mult * mju * (normal.x * tangent.y + normal.y * tangent.x);
      Scalar sigmayy = -mult * mju * (2 * normal.y * tangent.y);
      res.SetTension(sigmaxx, sigmayy, sigmaxy);
    } else
    {
      Scalar pSpeed = sqrt((lambda + Scalar(2.0) * mju) * invRho);

      Scalar omega = 2 * pi * pSpeed * invWaveLength;
      Scalar mult = cos(point * normal * k - omega * 0 // time 
                    + initialPhase);
      res.SetVelocity(normal * pSpeed * mult);

      Scalar sigmaxx = -mult * (lambda + 2 * mju * normal.x * normal.x);
      Scalar sigmaxy = -mult *           2 * mju * normal.x * normal.y ;
      Scalar sigmayy = -mult * (lambda + 2 * mju * normal.y * normal.y);
      res.SetTension(sigmaxx, sigmayy, sigmaxy);
    } */
    return res;
  }

private:
  Vector normal;
  Vector tangent;
  Scalar invWaveLength;
  Scalar initialPhase;
  bool shear;

  Scalar k;
}; 

template<typename ElasticSpace>
class ConstantPlaneWaveIniStateMaker: public IniStateMaker<ElasticSpace>
{
public:
  typedef typename ElasticSpace::Elastic                   Elastic;
  typedef typename ElasticSpace::SpaceType                 Space;
  typedef typename ElasticSpace::MediumParametersType      MediumParameters;
  SPACE_TYPEDEFS

  ConstantPlaneWaveIniStateMaker(const Vector& velocity, const Vector& center, Scalar waveLength, bool shear = false):
    normal(velocity.GetNorm()),
    center(center),
    waveLength(waveLength),
    velocity(velocity),
    shear(shear)
  {
    tangent = normal.GetPerpendicular();
  }

  void MakeParamsDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult) override
  {
    center     /= velocityDimensionlessMult;
    waveLength /= velocityDimensionlessMult;
    velocity   /= velocityDimensionlessMult;
  }

  Elastic GetValue(const Vector& point, const Scalar lambda, const Scalar mju, const Scalar invRho) override
  {
    Elastic res;

    Vector projection;
    GetProjection(point, center, normal, normal, projection);

    if ((point - projection).Len() > waveLength * Scalar(0.5))
    {
      res.SetZeroValues();
    } else
    {
      SetValues(Overload<Space>(), lambda, mju, invRho, res);
    }
    return res;
  }

  void SetValues(Overload<Space2>, Scalar lambda, Scalar mju, Scalar invRho, Elastic& res)
  {
    if (shear)
    {
      Scalar sSpeed = sqrt(mju * invRho);
      res.SetVelocity(tangent * velocity.Len());
      Scalar mult = velocity.Len() / sSpeed;

      Scalar sigmaxx = -mult * mju * (2 * normal.x * tangent.x);
      Scalar sigmaxy = -mult * mju * (normal.x * tangent.y + normal.y * tangent.x);
      Scalar sigmayy = -mult * mju * (2 * normal.y * tangent.y);
      res.SetTension(sigmaxx, sigmayy, sigmaxy);
      res.SetTension(0, 0, 0);
    } else
    {
      Scalar pSpeed = sqrt((lambda + Scalar(2.0) * mju) * invRho);

      res.SetVelocity(normal * velocity.Len());

      if (velocity.x > 0)
      {
        res.SetVelocity(-Vector::yAxis() * 0.05);
      } else
        res.SetVelocity(Vector::yAxis() * 0.05);

      Scalar mult = velocity.Len() / pSpeed;

      Scalar sigmaxx = -mult * (lambda + 2 * mju * normal.x * normal.x);
      Scalar sigmaxy = -mult * 2 * mju * normal.x * normal.y;
      Scalar sigmayy = -mult * (lambda + 2 * mju * normal.y * normal.y);
      res.SetTension(0, 0, 0);
    }
  }

  void SetValues(Overload<Space3>, Scalar lambda, Scalar mju, Scalar invRho, Elastic& res)
  {
    // TODO
    assert(0);
  }

  void GetProjection(const Space2::Vector& point, 
    const Space2::Vector& linePoint, const Space2::Vector& lineNormal, const Space2::Vector& projectionDirection,
    Space2::Vector& projectedPoint) const
  {
    ProjectPointAgainstLine(point, linePoint, lineNormal, projectionDirection, projectedPoint);
  }

  void GetProjection(const Space3::Vector& point,
    const Space3::Vector& planePoint, const Space3::Vector& planeNormal, const Space3::Vector& projectionDirection,
    Space3::Vector& projectedPoint) const
  {
    ProjectPointAgainstPlane(point, planePoint, planeNormal, projectionDirection, projectedPoint);
  }

private:
  Vector normal;
  Vector center;
  Vector tangent;
  Scalar waveLength;
  Vector velocity;
  Scalar initialPhase;
  bool shear;
}; 


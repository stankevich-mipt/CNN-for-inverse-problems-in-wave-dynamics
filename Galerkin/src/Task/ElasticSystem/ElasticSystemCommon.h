#pragma once

#include "../../Maths/Spaces.h"
#include "BoundaryInfos.h"
#include "ContactInfos.h"
#include "PointSources.h"

#include "SourceFunctors.h"
#include <vector>

#include "Eigen/Dense"

template <typename Space>
struct ElasticSystemBase;

template <>
struct ElasticSystemBase<Space2>
{
  static const int dimsCount = 5;
};

template <>
struct ElasticSystemBase<Space3>
{
  static const int dimsCount = 9;
};

template <typename Space>
struct ElasticSystemCommon: public ElasticSystemBase<Space>
{
  SPACE_TYPEDEFS
  typedef SourceFunctor<Space> SourceFunctorT;
  using ElasticSystemBase<Space>::dimsCount;

  #ifdef USE_DYNAMIC_MATRICIES
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXDim;
  #else
   typedef typename Eigen::Matrix<Scalar, dimsCount, dimsCount> MatrixXDim;
  #endif
  

  struct ValueTypeCommon
  {
    Scalar values[ElasticSystemBase<Space>::dimsCount];

    ValueTypeCommon() = default;
    virtual ~ValueTypeCommon() = default;

    ValueTypeCommon(Scalar initialValue)
    {
      std::fill_n(values, ElasticSystemBase<Space>::dimsCount, initialValue);
    }

    ValueTypeCommon& operator=(const ValueTypeCommon& other)
    {
      std::copy(other.values, other.values + ElasticSystemBase<Space>::dimsCount, values);
      return *this;
    }

    ValueTypeCommon& operator+=(const ValueTypeCommon& other)
    {
      for (IndexType i = 0; i < ElasticSystemBase<Space>::dimsCount; ++i)
        values[i] += other.values[i];
      return *this;
    }

    ValueTypeCommon& operator/=(Scalar value)
    {
      for (IndexType i = 0; i < ElasticSystemBase<Space>::dimsCount; ++i)
        values[i] /= value;
      return *this;
    }

    void SetZeroValues();
    Vector GetVelocity() const;
    void SetVelocity(const Vector& velocity);
    void SetTension(const Tensor& tension);

    void MakeDimension(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult);

    virtual Vector GetForce(const Vector& normal) const = 0;

    Scalar GetXX() const;
    Scalar GetYY() const;
    Scalar GetXY() const;
  };

  struct MediumParameters
  {
    SPACE_TYPEDEFS

    MediumParameters();
    MediumParameters(Scalar lambda, Scalar mju, Scalar invRho);

    void SetFromVelocities(Scalar pSpeed, Scalar sSpeed, Scalar rho);
    void SetFromLameParameters(Scalar lambda, Scalar mju, Scalar rho);
    void SetFromYoungModuleAndNju(Scalar youngModule, Scalar nju, Scalar rho);

    Scalar GetPSpeed() const;
    Scalar GetSSpeed() const;

    Scalar GetZp() const;
    Scalar GetZs() const;

    void MakeDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult);

    static const int ParamsCount = 3;
    union 
    {
      struct 
      {
        Scalar lambda;
        Scalar mju;
        Scalar invRho;
      };
      Scalar params[ParamsCount];
    };

    Scalar initialVolume;
    Scalar initialRho;

    bool IsZero() const;

    Vector flowVelocity;

    bool destroyed;
    bool dimensionless;

    struct Plasticity
    {
      // plasticity: 
      // s:s < 2 * k^2 
      // k = k0 + a * pressure
      Plasticity(): k0(std::numeric_limits<Scalar>::infinity()), a(0) /* without plasticity */, brittle(false),
        maxPlasticDeform(std::numeric_limits<Scalar>::infinity()),
        powderShearMult(0)
      {}
      Scalar k0;
      Scalar a;
      bool brittle;
      Scalar maxPlasticDeform;
      Scalar powderShearMult;
    };
    Plasticity plasticity;
    bool fixed;
  };
  
  void BuildXMatrix(const MediumParameters& mediumParameters, MatrixXDim& xMatrix);
  
  void BuildYMatrix(const MediumParameters& mediumParameters, MatrixXDim& yMatrix);
  
  void BuildRMatrix(const MediumParameters& interiorMediumParameters, 
                    const MediumParameters& exteriorMediumParameters, 
                    MatrixXDim& rMatrix);

  void BuildXnAuxMatrix(const MediumParameters& interiorMediumParameters, 
                        const MediumParameters& exteriorMediumParameters,
                        MatrixXDim& xnAuxMatrix);

  void BuildXnInteriorMatrix(const MediumParameters& interiorMediumParameters, 
                             const MediumParameters& exteriorMediumParameters,
                             const Vector& edgeNormal,
                             const MatrixXDim& xnAuxMatrix,
                             MatrixXDim& xnInteriorMatrix);

  void BuildXnExteriorMatrix(const MediumParameters& interiorMediumParameters, 
                             const MediumParameters& exteriorMediumParameters,
                             const Vector& edgeNormal, 
                             const MatrixXDim& xnAuxMatrix,
                             MatrixXDim& xnExteriorMatrix);

  void BuildBoundaryMatrix(IndexType interactionType, Eigen::Matrix<Scalar, 1, dimsCount>& boundaryMatrix);

  void BuildContactMatrices(IndexType interactionType,
    Eigen::Matrix<Scalar, 1, dimsCount>& leftContactMatrix,
    Eigen::Matrix<Scalar, 1, dimsCount>& rightContactMatrix);
  
  Scalar GetMaxWaveSpeed(const MediumParameters& mediumParameters) const;

  bool IsProperContact(const ValueTypeCommon& value0, const ValueTypeCommon& value1, 
    const Vector& contactNormal);

  BoundaryInfoFunctor<Space>* GetBoundaryInfoFunctor(IndexType interactionType);

  void SetFreeBoundary(IndexType interactionType, Scalar tensionDimensionlessMult, Scalar reflectionCoeff,
    VectorFunctor<Space>* externalForceFunctor = nullptr, IndexType dynamicContactInteractionType = IndexType(-1));

  void SetFixedBoundary(IndexType interactionType, Scalar velocityDimensionlessMult, Scalar reflectionCoeff,
    VectorFunctor<Space>* externalVelocityFunctor = nullptr);

  void SetAbsorbBoundary(IndexType interactiontype);
  void SetSymmetryBoundary(IndexType interactionType);
  void SetAntiSymmetryBoundary(IndexType interactionType);
  BoundaryConditions::Types GetBoundaryType(IndexType interactionType) const;

  void SetGlueContact(IndexType interactionType);
  void SetGlueContact(IndexType interactionType, Scalar maxShearStress, Scalar maxLongitudinalStress,
    IndexType dynamicBoundaryInteractionType);
  void SetGlideContact(IndexType interactionType);
  void SetFrictionContact(IndexType interactionType, Scalar frictionCoeff, IndexType dynamicBoundaryInteractionType);
  void SetFrictionContact(IndexType interactionType);
  void GetFrictionContactInfo(IndexType interactionType, Scalar& frictionCoeff);

  ContactConditions::Types GetContactType(IndexType interactionType) const;

  IndexType GetContactDynamicBoundaryType(IndexType contactInteractionType ) const;
  IndexType GetBoundaryDynamicContactType(IndexType boundaryInteractionType) const;

  void GetContactCriticalInfo(IndexType interactionType, 
    Scalar& maxShearStress, Scalar& maxLongitudinalStress);

  SourceFunctorT* GetSourceFunctor() const;
  void SetSourceFunctor(SourceFunctorT* sourceFunctor);
  void SetPointSources(const std::vector<PointSource<Space>* >& pointSources);
  std::vector< PointSource<Space>* > pointSources;

  class ZeroBoundaryInfoFunctor: public BoundaryInfoFunctor<Space> //default functor
  {
  public:
    SPACE_TYPEDEFS
    void operator()(const Vector&, const Vector&, const Scalar, Scalar* values) override
    {
      std::fill_n(values, dimsCount, Scalar(0.0));
    }
  };

  class FreeBoundaryInfoFunctor: public BoundaryInfoFunctor<Space>
  {
    SPACE_TYPEDEFS
  public:
    FreeBoundaryInfoFunctor(VectorFunctor<Space>* forceFunctor, Scalar tensionDimensionlessMult)
    {
      this->forceFunctor             = forceFunctor;
      this->tensionDimensionlessMult = tensionDimensionlessMult;
    }

    void operator()(const Vector& globalPoint, const Vector& externalNormal, const Scalar time, Scalar* values) override;

    void SetCurrentVelocity(const Vector& velocity) override
    {
      forceFunctor->SetCurrentVelocity(velocity);
    }

  private:
    VectorFunctor<Space>* forceFunctor;
    Scalar tensionDimensionlessMult; 
  };

  class FixedBoundaryInfoFunctor: public BoundaryInfoFunctor<Space>
  {
    SPACE_TYPEDEFS
  public:
    FixedBoundaryInfoFunctor(VectorFunctor<Space>* velocityFunctor, Scalar velocityDimensionlessMult)
    {
      this->velocityFunctor           = velocityFunctor;
      this->velocityDimensionlessMult = velocityDimensionlessMult;
    }

    void operator()(const Vector& globalPoint, const Vector& externalNormal, const Scalar time, Scalar* values) override
    {
      std::fill(values, values + dimsCount, 0);
      Vector externalVelocity = (*velocityFunctor)(globalPoint, externalNormal, time) / velocityDimensionlessMult;
      SetVelocity(externalVelocity, values);
    }

  private:
    void SetVelocity(const Vector& externalVelocity, Scalar* values);

    VectorFunctor<Space>* velocityFunctor;
    Scalar velocityDimensionlessMult;
  };

protected:
  void SetBoundaryDescription(IndexType interactionType, BoundaryDescription boundaryDescription);
  void SetContactDescription(IndexType interactionType, ContactDescription contactDescription);

  std::vector< FreeBoundaryInfo<Space>  > freeBoundaryInfos;
  std::vector< FixedBoundaryInfo<Space> > fixedBoundaryInfos;
  std::vector< BoundaryDescription >      boundaryDescriptions;

  std::vector< GlueContactInfo<Space> >     glueContactInfos;
  std::vector< FrictionContactInfo<Space> > frictionContactInfos;
  std::vector< ContactDescription >         contactDescriptions;
  SourceFunctorT* sourceFunctor; 

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include "ElasticSystemCommon.inl"

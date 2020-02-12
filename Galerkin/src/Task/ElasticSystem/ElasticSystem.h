#pragma once
#include "VectorFunctors.h"
#include "../../Maths/Spaces.h"
#include "ElasticSystemCommon.h"

template <typename Space>
struct ElasticSystem;

template <>
struct ElasticSystem<Space2>: public ElasticSystemCommon<Space2>
{
  SPACE2_TYPEDEFS
  typedef Space2                                        Space;
  typedef ElasticSystemCommon<Space>::ValueTypeCommon  ValueTypeCommon;
  typedef ElasticSystemCommon<Space>::MediumParameters MediumParameters;
  typedef ElasticSystemCommon<Space>::MatrixXDim MatrixXDim;

  using ElasticSystemCommon<Space2>::GetContactDynamicBoundaryType;

  struct ValueType: public ValueTypeCommon
  {
    using ValueTypeCommon::SetTension;
    using ValueTypeCommon::values;

    ValueType() : ValueTypeCommon()
    {
    }

    ValueType(Scalar initialValue) : ValueTypeCommon(initialValue)
    {
    }

    virtual ~ValueType() = default;

    void SetTension(Scalar xx, Scalar yy, Scalar xy);
    Scalar GetPressure() const;
    Scalar GetDeviatorSquare() const;
    Vector GetForce(const Vector& normal) const override;
    Tensor GetTension() const;
  };

  struct ElasticSpace: public Space2
  {
    SPACE2_TYPEDEFS
    typedef          Space2               SpaceType;
    typedef          ValueType            Elastic;
    typedef          MediumParameters     MediumParametersType;
  };

  ValueType GetRiemannSolution(const ValueType& interiorSolution, ValueType& exteriorSolution, 
    const MediumParameters& interiorParams, MediumParameters& exteriorParams,
    IndexType boundaryType, IndexType dynamicContactType);

  ValueType GetGlueRiemannSolution(const ValueType& interiorSolution, ValueType& exteriorSolution,
    const MediumParameters& interiorParams, const MediumParameters& exteriorParams);

  void BuildEdgeTransformMatrix(Vector edgeVertices[2], MatrixXDim& transformMatrix);
  void BuildEdgeTransformMatrixInv(Vector edgeVertices[2], MatrixXDim& transformMatrixInv);
};

template <>
struct ElasticSystem<Space3>: public ElasticSystemCommon<Space3>
{
  SPACE3_TYPEDEFS
  typedef Space3            Space;
  using ElasticSystemCommon<Space>::ValueTypeCommon;
  using ElasticSystemCommon<Space>::MediumParameters;
  using ElasticSystemCommon<Space>::GetContactDynamicBoundaryType;
  typedef ElasticSystemCommon<Space>::MatrixXDim MatrixXDim;

  struct ValueType: public ValueTypeCommon
  {
    ValueType() : ValueTypeCommon()
    {
    }

    ValueType(Scalar initialValue) : ValueTypeCommon(initialValue)
    {
    }

    using ValueTypeCommon::SetTension;
    using ValueTypeCommon::values;
    void SetTension(Scalar xx, Scalar yy, Scalar zz, Scalar xy, Scalar yz, Scalar xz);
    Tensor GetTension() const;

    Scalar GetZZ() const;
    Scalar GetYZ() const;
    Scalar GetXZ() const;
    Scalar GetPressure() const;
    Scalar GetDeviatorSquare() const;
    Vector GetForce(const Vector& normal) const override;
  };

  struct ElasticSpace: public Space3
  {
    SPACE3_TYPEDEFS
    typedef          Space3               SpaceType;
    typedef          ValueType            Elastic;
    typedef          MediumParameters     MediumParametersType;
  };

  ValueType GetRiemannSolution(const ValueType& interiorSolution, ValueType& exteriorSolution,
    const MediumParameters& interiorParams, MediumParameters& exteriorParams,
    IndexType boundaryType, IndexType dynamicContactType);

  ValueType GetGlueRiemannSolution(const ValueType& interiorSolution, ValueType& exteriorSolution,
    const MediumParameters& interiorParams, const MediumParameters& exteriorParams);

  void BuildFaceTransformMatrix(Vector faceVertices[3], MatrixXDim& transformMatrix);
  void BuildFaceTransformMatrixInv(Vector faceVertices[3], MatrixXDim& transformMatrixInv);
  void BuildZMatrix(const MediumParameters& mediumParameters, MatrixXDim& zMatrix);
};

#include "ElasticSystem.inl"

#include "../VolumeMethod/FunctionGetters/BoundaryFunctionGetter.h"

template <typename Space>
void ElasticSystemCommon<Space>::ValueTypeCommon::SetZeroValues()
{
  std::fill_n(values, dimsCount, Scalar(0.0));
}

template <>
Space2::Vector ElasticSystemCommon<Space2>::ValueTypeCommon::GetVelocity() const
{
  return Space2::Vector(values[3], values[4]);
}

template <>
Space3::Vector ElasticSystemCommon<Space3>::ValueTypeCommon::GetVelocity() const
{
  return Space3::Vector(values[6], values[7], values[8]);
}

template <>
void ElasticSystemCommon<Space2>::ValueTypeCommon::SetVelocity(const Vector& velocity)
{
  values[3] = velocity.x;
  values[4] = velocity.y;
}

template <>
void ElasticSystemCommon<Space3>::ValueTypeCommon::SetVelocity(const Vector& velocity)
{
  values[6] = velocity.x;
  values[7] = velocity.y;
  values[8] = velocity.z;
}

template <>
void ElasticSystemCommon<Space2>::ValueTypeCommon::SetTension(const Tensor& tension)
{
  values[0] = tension.xx;
  values[1] = tension.yy;
  values[2] = tension.xy;
}

template <>
void ElasticSystemCommon<Space3>::ValueTypeCommon::SetTension(const Tensor& tension)
{
  values[0] = tension.xx;
  values[1] = tension.yy;
  values[2] = tension.zz;
  values[3] = tension.xy;
  values[4] = tension.yz;
  values[5] = tension.xz;
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::ValueTypeCommon::GetXX() const
{
  return values[0];
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::ValueTypeCommon::GetYY() const
{
  return values[1];
}

template <>
Space2::Scalar ElasticSystemCommon<Space2>::ValueTypeCommon::GetXY() const
{
  return values[2];
}

template <>
Space3::Scalar ElasticSystemCommon<Space3>::ValueTypeCommon::GetXY() const
{
  return values[3];
}

template <>
void ElasticSystemCommon<Space2>::ValueTypeCommon::MakeDimension(
  Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult)
{
  for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
  {
    switch (valueIndex)
    {
      case 0: case 1: case 2:
        values[valueIndex] *= tensionDimensionlessMult;
      break;
      case 3: case 4:
        values[valueIndex] *= velocityDimensionlessMult;
      break;
      default:
        assert(0);
      break;
    }
  }
}

template <>
void ElasticSystemCommon<Space3>::ValueTypeCommon::MakeDimension(
  Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult)
{
  for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
  {
    switch (valueIndex)
    {
      case 0: case 1: case 2: case 3: case 4: case 5:
        values[valueIndex] *= tensionDimensionlessMult;
      break;
      case 6: case 7: case 8:
        values[valueIndex] *= velocityDimensionlessMult;
      break;
      default:
        assert(0);
      break;
    }
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildXMatrix(const MediumParameters& mediumParameters, MatrixXDim& xMatrix)
{
  xMatrix << 
    mediumParameters.flowVelocity.x,                               0,                               0, -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju),                               0,
                                  0, mediumParameters.flowVelocity.x,                               0,                                        -mediumParameters.lambda,                               0,
                                  0,                               0, mediumParameters.flowVelocity.x,                                                               0,           -mediumParameters.mju,
           -mediumParameters.invRho,                               0,                               0,                                 mediumParameters.flowVelocity.x,                               0,
                                  0,                               0,        -mediumParameters.invRho,                                                               0, mediumParameters.flowVelocity.x;
}

template <>
void ElasticSystemCommon<Space3>::BuildXMatrix(const MediumParameters& mediumParameters, MatrixXDim& xMatrix)
{
  xMatrix <<
    mediumParameters.flowVelocity.x,                               0,                               0,                               0,                               0,                               0, -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju),                               0,                               0,
                                  0, mediumParameters.flowVelocity.x,                               0,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,                               0,                               0,
                                  0,                               0, mediumParameters.flowVelocity.x,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,                               0,                               0,
                                  0,                               0,                               0, mediumParameters.flowVelocity.x,                               0,                               0,                                                               0,           -mediumParameters.mju,                               0,
                                  0,                               0,                               0,                               0, mediumParameters.flowVelocity.x,                               0,                                                               0,                               0,                               0,
                                  0,                               0,                               0,                               0,                               0, mediumParameters.flowVelocity.x,                                                               0,                               0,           -mediumParameters.mju,
           -mediumParameters.invRho,                               0,                               0,                               0,                               0,                               0,                                 mediumParameters.flowVelocity.x,                               0,                               0,
                                  0,                               0,                               0,        -mediumParameters.invRho,                               0,                               0,                                                               0, mediumParameters.flowVelocity.x,                               0,
                                  0,                               0,                               0,                               0,                               0,        -mediumParameters.invRho,                                                               0,                               0, mediumParameters.flowVelocity.x;
}

template <>
void ElasticSystemCommon<Space2>::BuildYMatrix(const MediumParameters& mediumParameters, MatrixXDim& yMatrix)
{
  yMatrix <<
    mediumParameters.flowVelocity.y,                               0,                               0,                               0,                                        -mediumParameters.lambda,
                                  0, mediumParameters.flowVelocity.y,                               0,                               0, -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju),
                                  0,                               0, mediumParameters.flowVelocity.y,           -mediumParameters.mju,                                                               0,
                                  0,                               0,        -mediumParameters.invRho, mediumParameters.flowVelocity.y,                                                               0,
                                  0,        -mediumParameters.invRho,                               0,                               0,                                 mediumParameters.flowVelocity.y;
}

template <>
void ElasticSystemCommon<Space3>::BuildYMatrix(const MediumParameters& mediumParameters, MatrixXDim& yMatrix)
{
  yMatrix <<
    mediumParameters.flowVelocity.y,                               0,                               0,                               0,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,                               0,
                                  0, mediumParameters.flowVelocity.y,                               0,                               0,                               0,                               0,                               0, -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju),                               0,
                                  0,                               0, mediumParameters.flowVelocity.y,                               0,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,                               0,
                                  0,                               0,                               0, mediumParameters.flowVelocity.y,                               0,                               0,           -mediumParameters.mju,                                                               0,                               0,
                                  0,                               0,                               0,                               0, mediumParameters.flowVelocity.y,                               0,                               0,                                                               0,           -mediumParameters.mju,
                                  0,                               0,                               0,                               0,                               0, mediumParameters.flowVelocity.y,                               0,                                                               0,                               0,
                                  0,                               0,                               0,        -mediumParameters.invRho,                               0,                               0, mediumParameters.flowVelocity.y,                                                               0,                               0,
                                  0,        -mediumParameters.invRho,                               0,                               0,                               0,                               0,                               0,                                 mediumParameters.flowVelocity.y,                               0,
                                  0,                               0,                               0,                               0,        -mediumParameters.invRho,                               0,                               0,                                                               0, mediumParameters.flowVelocity.y;
}

template <>
void ElasticSystemCommon<Space2>::BuildRMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters, 
  MatrixXDim& rMatrix)
{
  rMatrix << 
    interiorMediumParameters.lambda + 2 * interiorMediumParameters.mju,                                    0, 0,                                     0, exteriorMediumParameters.lambda + 2 * exteriorMediumParameters.mju,
    interiorMediumParameters.lambda                                   ,                                    0, 1,                                     0, exteriorMediumParameters.lambda                                   ,
                                                                     0,         interiorMediumParameters.mju, 0,          exteriorMediumParameters.mju,                                                                  0,
                                  interiorMediumParameters.GetPSpeed(),                                    0, 0,                                     0,                              -exteriorMediumParameters.GetPSpeed(),
                                                                     0, interiorMediumParameters.GetSSpeed(), 0, -exteriorMediumParameters.GetSSpeed(),                                                                  0;
}

template <>
void ElasticSystemCommon<Space3>::BuildRMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters, 
  MatrixXDim& rMatrix)
{
  rMatrix <<
    interiorMediumParameters.lambda + 2 * interiorMediumParameters.mju,                                    0,                                    0, 0, 0, 0,                                     0,                                     0, exteriorMediumParameters.lambda + 2 * exteriorMediumParameters.mju,
    interiorMediumParameters.lambda                                   ,                                    0,                                    0, 0, 1, 0,                                     0,                                     0, exteriorMediumParameters.lambda                                   ,
    interiorMediumParameters.lambda                                   ,                                    0,                                    0, 0, 0, 1,                                     0,                                     0, exteriorMediumParameters.lambda                                   ,
                                                                     0,         interiorMediumParameters.mju,                                    0, 0, 0, 0,                                     0,          exteriorMediumParameters.mju,                                                                  0,
                                                                     0,                                    0,                                    0, 1, 0, 0,                                     0,                                     0,                                                                  0,
                                                                     0,                                    0,         interiorMediumParameters.mju, 0, 0, 0,          exteriorMediumParameters.mju,                                     0,                                                                  0,
                                  interiorMediumParameters.GetPSpeed(),                                    0,                                    0, 0, 0, 0,                                     0,                                     0,                              -exteriorMediumParameters.GetPSpeed(),
                                                                     0, interiorMediumParameters.GetSSpeed(),                                    0, 0, 0, 0,                                     0, -exteriorMediumParameters.GetSSpeed(),                                                                  0,
                                                                     0,                                    0, interiorMediumParameters.GetSSpeed(), 0, 0, 0, -exteriorMediumParameters.GetSSpeed(),                                     0,                                                                  0;
}

template <>
void ElasticSystemCommon<Space2>::BuildXnAuxMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  MatrixXDim& xnAuxMatrix)
{  
/*
  Eigen::Matrix<Scalar, dimsCount, dimsCount> rMatrix;
  BuildRMatrix(interiorMediumParameters, exteriorMediumParameters, rMatrix);

  Eigen::Matrix<Scalar, 1, dimsCount> lambdaMatrix;
  lambdaMatrix << 
    interiorMediumParameters.GetPSpeed(), 
    interiorMediumParameters.GetSSpeed(), 0, 0, 0;
  
  xnAuxMatrix.noalias() = rMatrix * lambdaMatrix.asDiagonal() * rMatrix.inverse();
*/

  Scalar cpInt = interiorMediumParameters.GetPSpeed();
  Scalar cpExt = exteriorMediumParameters.GetPSpeed();
  Scalar csInt = interiorMediumParameters.GetSSpeed();
  Scalar csExt = exteriorMediumParameters.GetSSpeed();

  Scalar rhoInt = 1 / interiorMediumParameters.invRho;
  Scalar rhoExt = 1 / exteriorMediumParameters.invRho;

  Scalar zpInt = cpInt * rhoInt;
  Scalar zsInt = csInt * rhoInt;
  Scalar zpExt = cpExt * rhoExt;
  Scalar zsExt = csExt * rhoExt;

  Scalar coeffP = 1 / (zpInt + zpExt);
  Scalar coeffS = Scalar(1.0) / (rhoInt + rhoExt);
  
  if (zsInt + zsExt > std::numeric_limits<Scalar>::epsilon())
  {
    coeffS = csInt / (zsInt + zsExt);
  }

  xnAuxMatrix << 
                      coeffP * zpInt * cpInt, 0,              0,                   coeffP * zpInt * zpExt * cpInt,                      0,
    coeffP * interiorMediumParameters.lambda, 0,              0, coeffP * interiorMediumParameters.lambda * zpExt,                      0,
                                           0, 0, coeffS * zsInt,                                                0, coeffS * zsInt * zsExt,
                              coeffP * cpInt, 0,              0,                           coeffP * cpInt * zpExt,                      0,
                                           0, 0,         coeffS,                                                0,         coeffS * zsExt;
}

template <>
void ElasticSystemCommon<Space3>::BuildXnAuxMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  MatrixXDim& xnAuxMatrix)
{
/*
  Eigen::Matrix<Scalar, dimsCount, dimsCount> rMatrix;
  BuildRMatrix(interiorMediumParameters, exteriorMediumParameters, rMatrix);

  Eigen::Matrix<Scalar, 1, dimsCount> lambdaMatrix;
  lambdaMatrix << interiorMediumParameters.GetPSpeed(),
                  interiorMediumParameters.GetSSpeed(),
                  interiorMediumParameters.GetSSpeed(),
                  0, 0, 0, 0, 0, 0;

  xnAuxMatrix.noalias() = rMatrix * lambdaMatrix.asDiagonal() * rMatrix.inverse();
*/

  Scalar cpInt = interiorMediumParameters.GetPSpeed();
  Scalar cpExt = exteriorMediumParameters.GetPSpeed();
  Scalar csInt = interiorMediumParameters.GetSSpeed();
  Scalar csExt = exteriorMediumParameters.GetSSpeed();

  Scalar rhoInt = 1 / interiorMediumParameters.invRho;
  Scalar rhoExt = 1 / exteriorMediumParameters.invRho;

  Scalar zpInt = cpInt * rhoInt;
  Scalar zsInt = csInt * rhoInt;
  Scalar zpExt = cpExt * rhoExt;
  Scalar zsExt = csExt * rhoExt;

  Scalar coeffP = 1 / (zpInt + zpExt);
  Scalar coeffS = Scalar(1.0) / (rhoInt + rhoExt);

  if (zsInt + zsExt > std::numeric_limits<Scalar>::epsilon())
  {
    coeffS = csInt / (zsInt + zsExt);
  }

  xnAuxMatrix <<
                      coeffP * zpInt * cpInt, 0, 0,              0, 0,              0,                   coeffP * zpInt * zpExt * cpInt,                      0,                      0,
    coeffP * interiorMediumParameters.lambda, 0, 0,              0, 0,              0, coeffP * zpExt * interiorMediumParameters.lambda,                      0,                      0,
    coeffP * interiorMediumParameters.lambda, 0, 0,              0, 0,              0, coeffP * zpExt * interiorMediumParameters.lambda,                      0,                      0,
                                           0, 0, 0, coeffS * zsInt, 0,              0,                                                0, coeffS * zsInt * zsExt,                      0,
                                           0, 0, 0,              0, 0,              0,                                                0,                      0,                      0,
                                           0, 0, 0,              0, 0, coeffS * zsInt,                                                0,                      0, coeffS * zsInt * zsExt,
                              coeffP * cpInt, 0, 0,              0, 0,              0,                           coeffP * zpExt * cpInt,                      0,                      0,
                                           0, 0, 0,         coeffS, 0,              0,                                                0,         coeffS * zsExt,                      0,
                                           0, 0, 0,              0, 0,         coeffS,                                                0,                      0,         coeffS * zsExt;
}

template <>
void ElasticSystemCommon<Space2>::BuildXnInteriorMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  const Vector& edgeNormal, 
  const MatrixXDim& xnAuxMatrix,
  MatrixXDim& xnInteriorMatrix)
{
  // BuildXnAuxMatrix(interiorMediumParameters, exteriorMediumParameters, xnInteriorMatrix);
  xnInteriorMatrix = xnAuxMatrix;

  // XMatrix
  xnInteriorMatrix(0, 3) -= interiorMediumParameters.lambda + Scalar(2.0) * interiorMediumParameters.mju;
  xnInteriorMatrix(1, 3) -= interiorMediumParameters.lambda;
  xnInteriorMatrix(2, 4) -= interiorMediumParameters.mju;
  xnInteriorMatrix(3, 0) -= interiorMediumParameters.invRho;
  xnInteriorMatrix(4, 2) -= interiorMediumParameters.invRho;

  Scalar un = interiorMediumParameters.flowVelocity * edgeNormal;
  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    xnInteriorMatrix(row, row) += (un + fabs(un)) * Scalar(0.5);
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildXnInteriorMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  const Vector& faceNormal, 
  const MatrixXDim& xnAuxMatrix,
  MatrixXDim& xnInteriorMatrix)
{
  //BuildXnAuxMatrix(interiorMediumParameters, exteriorMediumParameters, xnInteriorMatrix);
  xnInteriorMatrix = xnAuxMatrix;

  // XMatrix
  xnInteriorMatrix(0, 6) -= interiorMediumParameters.lambda + Scalar(2.0) * interiorMediumParameters.mju;
  xnInteriorMatrix(1, 6) -= interiorMediumParameters.lambda;
  xnInteriorMatrix(2, 6) -= interiorMediumParameters.lambda;
  xnInteriorMatrix(3, 7) -= interiorMediumParameters.mju;
  xnInteriorMatrix(5, 8) -= interiorMediumParameters.mju;
  xnInteriorMatrix(6, 0) -= interiorMediumParameters.invRho;
  xnInteriorMatrix(7, 3) -= interiorMediumParameters.invRho;
  xnInteriorMatrix(8, 5) -= interiorMediumParameters.invRho;

  Scalar un = interiorMediumParameters.flowVelocity * faceNormal;
  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    xnInteriorMatrix(row, row) += (un + fabs(un)) * Scalar(0.5);
  }
}

template <typename Space>
void ElasticSystemCommon<Space>::BuildXnExteriorMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  const Vector& normal, 
  const MatrixXDim& xnAuxMatrix,
  MatrixXDim& xnExteriorMatrix)
{
  //BuildXnAuxMatrix(interiorMediumParameters, exteriorMediumParameters, xnExteriorMatrix);
  //xnExteriorMatrix *= Scalar(-1.0);
  xnExteriorMatrix = -xnAuxMatrix;

  Scalar un = interiorMediumParameters.flowVelocity * normal;
  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    xnExteriorMatrix(row, row) += (un - fabs(un)) * Scalar(0.5);
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildBoundaryMatrix(IndexType interactionType, Eigen::Matrix<Scalar, 1, dimsCount>& boundaryMatrix)
{
  IndexType boundaryType = boundaryDescriptions[interactionType].type;
  Scalar reflectionCoeff = boundaryDescriptions[interactionType].reflectionCoeff;

  switch(boundaryType)
  {
    case BoundaryConditions::Absorb:
    {
      boundaryMatrix << 0, 0, 0, 0, 0;
    } break;
    case BoundaryConditions::Free:
    {
      boundaryMatrix << -reflectionCoeff, 
                         reflectionCoeff, 
                        -reflectionCoeff, 
                         reflectionCoeff, 
                         reflectionCoeff;
    } break;
    case BoundaryConditions::Fixed:
    {
      boundaryMatrix << reflectionCoeff, 
                        reflectionCoeff, 
                        reflectionCoeff, 
                       -reflectionCoeff, 
                       -reflectionCoeff;
    } break;
    case BoundaryConditions::Symmetry:
    {
      boundaryMatrix << Scalar(1.0),
                        Scalar(1.0), 
                       -Scalar(1.0), // sigma.xy 
                       -Scalar(1.0), // v.x 
                        Scalar(1.0);
    } break;
    case BoundaryConditions::AntiSymmetry:
    {
      boundaryMatrix << -Scalar(1.0), // sigma.xx 
                         Scalar(1.0), 
                         Scalar(1.0), 
                         Scalar(1.0), 
                        -Scalar(1.0); // v.y
    } break;
    default:
      assert(0);
    break;
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildBoundaryMatrix(IndexType interactionType, Eigen::Matrix<Scalar, 1, dimsCount>& boundaryMatrix)
{
  IndexType boundaryType = boundaryDescriptions[interactionType].type;
  Scalar reflectionCoeff = boundaryDescriptions[interactionType].reflectionCoeff;

  switch (boundaryType)
  {
    case BoundaryConditions::Absorb:
    {
      boundaryMatrix << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    } break;
    case BoundaryConditions::Free:
    {
      boundaryMatrix << 
        -reflectionCoeff, // sigma.xx
         reflectionCoeff,
         reflectionCoeff,
        -reflectionCoeff, // sigma.xy
         reflectionCoeff,
        -reflectionCoeff, // sigma.xz
         reflectionCoeff,
         reflectionCoeff,
         reflectionCoeff;
    } break;
    case BoundaryConditions::Fixed:
    {
      boundaryMatrix <<
        reflectionCoeff,
        reflectionCoeff,
        reflectionCoeff,
        reflectionCoeff,
        reflectionCoeff,
        reflectionCoeff,
       -reflectionCoeff, // v.x
       -reflectionCoeff, // v.y
       -reflectionCoeff; // v.z
    } break;
    case BoundaryConditions::Symmetry:
    {
      boundaryMatrix <<
        reflectionCoeff,
        reflectionCoeff,
        reflectionCoeff,
       -reflectionCoeff, // sigma xy
       -reflectionCoeff, // sigma yz
       -reflectionCoeff, // sigma xz
       -reflectionCoeff, // v.x
        reflectionCoeff,
        reflectionCoeff;
    } break;
    case BoundaryConditions::AntiSymmetry:
    {
      boundaryMatrix <<
        -reflectionCoeff, // sigma.xx
         reflectionCoeff,
         reflectionCoeff,
         reflectionCoeff,
         reflectionCoeff,
         reflectionCoeff,
         reflectionCoeff,
        -reflectionCoeff, // v.y
        -reflectionCoeff; // v.z
    } break;
    default:
      assert(0);
    break;
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildContactMatrices(IndexType interactionType, 
  Eigen::Matrix<Scalar, 1, dimsCount>& leftContactMatrix, Eigen::Matrix<Scalar, 1, dimsCount>& rightContactMatrix)
{
  IndexType contactType = contactDescriptions[interactionType].type;

  switch(contactType)
  {
    case ContactConditions::Glue:
    {
      leftContactMatrix  << 0, 0, 0, 0, 0;
      rightContactMatrix << 1, 1, 1, 1, 1;
    } break;
    case ContactConditions::Glide:
    {
      // it`s correct when materials on the both sides of contact are the same 
      leftContactMatrix  << 0, 1, -1, 0, 1;
      rightContactMatrix << 1, 0,  0, 1, 0;
    } break;
    case ContactConditions::Friction:
    {
      // friction should be computed as a dynamic contact
      assert(0);
    } break;
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildContactMatrices(IndexType interactionType,
  Eigen::Matrix<Scalar, 1, dimsCount>& leftContactMatrix, Eigen::Matrix<Scalar, 1, dimsCount>& rightContactMatrix)
{
  IndexType contactType = contactDescriptions[interactionType].type;

  switch (contactType)
  {
    case ContactConditions::Glue:
    {
      leftContactMatrix  << 0, 0, 0, 0, 0, 0, 0, 0, 0;
      rightContactMatrix << 1, 1, 1, 1, 1, 1, 1, 1, 1;
    } break;
    case ContactConditions::Glide:
    {
      // TODO
      leftContactMatrix  << 0, 1, 1, -1, -1, -1, 0, 1, 1;
      rightContactMatrix << 1, 0, 0,  0,  0,  0, 1, 0, 0;
    } break;
    case ContactConditions::Friction:
    {
      // friction should be computed as a dynamic contact
      assert(0);
    } break;
    default:
      assert(0);
    break;
  }
}

template<typename Space>
void ElasticSystemCommon<Space>::SetAbsorbBoundary(IndexType interactionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Absorb;
  newBoundary.infoIndex = -1;
  SetBoundaryDescription(interactionType, newBoundary);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetSymmetryBoundary(IndexType interactionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Symmetry;
  newBoundary.infoIndex = -1;
  SetBoundaryDescription(interactionType, newBoundary);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetAntiSymmetryBoundary(IndexType interactionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::AntiSymmetry;
  newBoundary.infoIndex = -1;
  SetBoundaryDescription(interactionType, newBoundary);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetGlueContact(IndexType interactionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Glue;
  newContact.infoIndex = -1;
  SetContactDescription(interactionType, newContact);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetGlueContact(IndexType interactionType, 
  Scalar maxShearStress, Scalar maxLongitudinalStress, IndexType dynamicBoundaryInteractionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Glue;
  newContact.infoIndex = glueContactInfos.size();
  SetContactDescription(interactionType, newContact);

  glueContactInfos.push_back(
    GlueContactInfo<Space>(maxShearStress, maxLongitudinalStress, dynamicBoundaryInteractionType));
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFrictionContact(IndexType interactionType, Scalar frictionCoeff, IndexType dynamicBoundaryInteractionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Friction;
  newContact.infoIndex = frictionContactInfos.size();
  SetContactDescription(interactionType, newContact);

  frictionContactInfos.push_back(FrictionContactInfo<Space>(frictionCoeff, dynamicBoundaryInteractionType));
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFrictionContact(IndexType interactionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Friction;
  newContact.infoIndex = frictionContactInfos.size();
  SetContactDescription(interactionType, newContact);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetGlideContact(IndexType interactionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Glide;
  newContact.infoIndex = -1;
  SetContactDescription(interactionType, newContact);
}

template<typename Space>
typename Space::IndexType ElasticSystemCommon<Space>::GetContactDynamicBoundaryType(IndexType contactInteractionType) const
{
  if(contactInteractionType < contactDescriptions.size() && 
     contactDescriptions[contactInteractionType].type == ContactConditions::Glue &&
     contactDescriptions[contactInteractionType].infoIndex != IndexType(-1))
  {
    IndexType glueInfoIndex = contactDescriptions[contactInteractionType].infoIndex;
    return glueInfoIndex != IndexType(-1) ? glueContactInfos.at(glueInfoIndex).dynamicBoundaryInteractionType : IndexType(-1);
  }
  return IndexType(-1);
}

template<typename Space>
typename Space::IndexType ElasticSystemCommon<Space>::GetBoundaryDynamicContactType(IndexType boundaryInteractionType) const
{
  if(boundaryInteractionType < boundaryDescriptions.size() &&
     boundaryDescriptions[boundaryInteractionType].type == BoundaryConditions::Free &&
     boundaryDescriptions[boundaryInteractionType].infoIndex != IndexType(-1))
  {
    IndexType freeInfoIndex = boundaryDescriptions[boundaryInteractionType].infoIndex;
    return freeInfoIndex != IndexType(-1) ? freeBoundaryInfos[freeInfoIndex].dynamicContactInteractionType : IndexType(-1);
  }
  return IndexType(-1);
}

template<typename Space>
void ElasticSystemCommon<Space>::GetFrictionContactInfo(IndexType interactionType, Scalar &frictionCoeff)
{
  if(contactDescriptions[interactionType].type == ContactConditions::Friction &&
     contactDescriptions[interactionType].infoIndex != IndexType(-1))
  {
    frictionCoeff = frictionContactInfos[contactDescriptions[interactionType].infoIndex].frictionCoeff;
    return;
  }
  assert(0);
  frictionCoeff = 0; //should never happen. hopefully.
}

template<typename Space>
void ElasticSystemCommon<Space>::GetContactCriticalInfo(IndexType interactionType, 
  Scalar &maxShearStress, Scalar &maxLongitudinalStress)
{
  assert(GetContactDynamicBoundaryType(interactionType) != IndexType(-1));
  if(contactDescriptions[interactionType].type == ContactConditions::Glue &&
     contactDescriptions[interactionType].infoIndex != IndexType(-1))
  {
    maxShearStress        = glueContactInfos[contactDescriptions[interactionType].infoIndex].maxShearStress;
    maxLongitudinalStress = glueContactInfos[contactDescriptions[interactionType].infoIndex].maxLongitudinalStress;
    return;
  }
  assert(0);
  maxShearStress = 0; //should never happen. hopefully.
  maxLongitudinalStress = 0;
}

template<typename Space>
BoundaryConditions::Types ElasticSystemCommon<Space>::GetBoundaryType(IndexType interactionType) const
{
  return boundaryDescriptions[interactionType].type;
}

template<typename Space>
ContactConditions::Types ElasticSystemCommon<Space>::GetContactType(IndexType interactionType) const
{
  return contactDescriptions[interactionType].type;
}

template<typename Space>
typename ElasticSystemCommon<Space>::SourceFunctorT* ElasticSystemCommon<Space>::GetSourceFunctor() const
{
  return sourceFunctor;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetSourceFunctor(SourceFunctorT* sourceFunctor)
{
  this->sourceFunctor = sourceFunctor;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetPointSources(const std::vector<PointSource<Space>* >& pointSources)
{
  this->pointSources = pointSources;
}

template<typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::GetMaxWaveSpeed(const MediumParameters& mediumParameters) const
{
  return mediumParameters.GetPSpeed();
}

template<typename Space>
void ElasticSystemCommon<Space>::SetBoundaryDescription(IndexType interactionType, BoundaryDescription boundaryDescription)
{
  if(boundaryDescriptions.size() < interactionType + 1)
  {
    boundaryDescriptions.resize(interactionType + 1);
  }
  boundaryDescriptions[interactionType] = boundaryDescription;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetContactDescription(IndexType interactionType, ContactDescription contactDescription)
{
  if(contactDescriptions.size() < interactionType + 1)
  {
    contactDescriptions.resize(interactionType + 1);
  }
  contactDescriptions[interactionType] = contactDescription;
}

template <typename Space>
ElasticSystemCommon<Space>::MediumParameters::MediumParameters(): 
  lambda(Scalar(2.0)), 
  mju(Scalar(1.0)), 
  invRho(Scalar(1.0)),
  destroyed(false),
  dimensionless(false),
  fixed(false)
{
}

template <typename Space>
ElasticSystemCommon<Space>::MediumParameters::MediumParameters(
  Scalar lambda, Scalar mju, Scalar invRho): 
  lambda(lambda), 
  mju(mju), 
  invRho(invRho),
  destroyed(false),
  dimensionless(false),
  fixed(false)
{
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::SetFromVelocities(Scalar pSpeed, Scalar sSpeed, Scalar rho)
{
  mju = Sqr(sSpeed) * rho;
  lambda = (Sqr(pSpeed) - 2.0 * Sqr(sSpeed)) * rho;
  invRho = 1 / rho;
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::SetFromLameParameters(Scalar lambda, Scalar mju, Scalar rho)
{
  this->lambda = lambda;
  this->mju    = mju;
  invRho = 1 / rho;
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::SetFromYoungModuleAndNju(Scalar youngModule, Scalar nju, Scalar rho)
{
  this->lambda = nju * youngModule / (1 + nju) / (1 - 2 * nju);
  this->mju    = youngModule / Scalar(2.0) / (1 + nju);
  invRho = 1 / rho;
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetPSpeed() const
{
  return sqrt((lambda + Scalar(2.0) * mju) * invRho);
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetSSpeed() const
{
  return sqrt(mju * invRho);
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetZp() const
{
  return GetPSpeed() / invRho;
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetZs() const
{
  return GetSSpeed() / invRho;
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::MakeDimensionless(
  Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult)
{
  if (!dimensionless)
  {
    dimensionless = true;
    lambda        /= tensionDimensionlessMult;
    mju           /= tensionDimensionlessMult;
    invRho        *= tensionDimensionlessMult / Sqr(velocityDimensionlessMult);
    plasticity.k0 /= tensionDimensionlessMult;
    flowVelocity  /= velocityDimensionlessMult;
  }
}

template <typename Space>
bool ElasticSystemCommon<Space>::MediumParameters::IsZero() const
{
  bool isZero = true;
  for (IndexType paramIndex = 0; paramIndex < ParamsCount; ++paramIndex)
  {
    if (params[paramIndex] > std::numeric_limits<Scalar>::epsilon()) isZero = false;
  }
  return isZero;
}

template<typename Space>
bool ElasticSystemCommon<Space>::IsProperContact(const ValueTypeCommon& value0, 
                                                 const ValueTypeCommon& value1, const Vector& contactNormal)
{
  // return true;
  Scalar eps = std::numeric_limits<Scalar>::epsilon();

  if (value1.GetForce(-contactNormal) * contactNormal - value0.GetForce(contactNormal) * contactNormal > -eps &&
    (value0.GetVelocity() - value1.GetVelocity()) * contactNormal > eps) return true;
  return false;
}

template<typename Space>
BoundaryInfoFunctor<Space>*
  ElasticSystemCommon<Space>::GetBoundaryInfoFunctor(IndexType interactionType)
{
  int infoIndex = boundaryDescriptions[interactionType].infoIndex;

  BoundaryFunctionGetter<Space>* getter = nullptr;

  switch(boundaryDescriptions[interactionType].type)
  {
    case BoundaryConditions::Absorb:
    {
      return nullptr;
    }break;
    case BoundaryConditions::Free:
    {
      assert(infoIndex != -1);
      return freeBoundaryInfos[infoIndex].boundaryFunctor;
    }break;
    case BoundaryConditions::Fixed:
    {
      assert(infoIndex != -1);
      return fixedBoundaryInfos[infoIndex].boundaryFunctor;
    }break;
    case BoundaryConditions::Symmetry:
    {
      return nullptr;
    }break;
    case BoundaryConditions::AntiSymmetry:
    {
      return nullptr;
    }break;
    default:
      assert(0);
    break;
  }
  assert(0); //unknown interation type
  return nullptr;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFreeBoundary(IndexType interactionType, Scalar tensionDimensionlessMult, Scalar reflectionCoeff,
  VectorFunctor<Space>* externalForceFunctor, 
  IndexType dynamicContactInteractionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Free;
  newBoundary.infoIndex = freeBoundaryInfos.size();
  newBoundary.reflectionCoeff = reflectionCoeff;
  SetBoundaryDescription(interactionType, newBoundary);

  FreeBoundaryInfo<Space> freeInfo;

  if (externalForceFunctor)
    freeInfo.boundaryFunctor = externalForceFunctor ? new FreeBoundaryInfoFunctor(externalForceFunctor, tensionDimensionlessMult) : nullptr;
  freeInfo.dynamicContactInteractionType = dynamicContactInteractionType;
  freeBoundaryInfos.push_back(freeInfo);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFixedBoundary(IndexType interactionType, Scalar velocityDimensionlessMult, Scalar reflectionCoeff,
  VectorFunctor<Space>* externalVelocityFunctor)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Fixed;
  newBoundary.infoIndex = fixedBoundaryInfos.size();
  newBoundary.reflectionCoeff = reflectionCoeff;
  SetBoundaryDescription(interactionType, newBoundary);
  
  FixedBoundaryInfo<Space> fixedInfo;
  fixedInfo.boundaryFunctor = externalVelocityFunctor ? new FixedBoundaryInfoFunctor(externalVelocityFunctor, velocityDimensionlessMult) : nullptr;

  fixedBoundaryInfos.push_back(fixedInfo);
}

template <>
void ElasticSystemCommon<Space2>::FreeBoundaryInfoFunctor::operator()(
  const Vector& globalPoint, 
  const Vector& externalNormal, 
  const Scalar time, Scalar* values)
{
  Vector externalForce = (*forceFunctor)(globalPoint, externalNormal, time) / tensionDimensionlessMult;
  Vector force(externalForce * externalNormal, externalNormal ^ externalForce);

  /*const Scalar frictionCoeff = Scalar(0.002);
  force.y *= -frictionCoeff * externalForce.GetNorm() * externalNormal;*/

  Tensor tension(force.x, force.y, 0);
  Tensor rotatedTension = tension.RotateAxes(-externalNormal);

  values[0] = rotatedTension.xx;
  values[1] = rotatedTension.yy;
  values[2] = rotatedTension.xy;
  values[3] = values[4] = 0;
}

template <>
void ElasticSystemCommon<Space3>::FreeBoundaryInfoFunctor::operator()(
  const Vector& globalPoint,
  const Vector& externalNormal,
  const Scalar time, Scalar* values)
{
  Vector externalForce = (*forceFunctor)(globalPoint, externalNormal, time) / tensionDimensionlessMult;

  Vector tangent0 = (externalNormal ^ externalForce).GetNorm();
  Vector tangent1 = (externalNormal ^ tangent0).GetNorm();

  Vector force(externalForce * externalNormal, externalForce * tangent0, externalForce * tangent1);

  Tensor tension(force.x, force.y, force.z, 0, 0, 0);
  Tensor rotatedTension = tension.RotateAxes(externalNormal, tangent0, tangent1);

  values[0] = rotatedTension.xx;
  values[1] = rotatedTension.yy;
  values[2] = rotatedTension.zz;
  values[3] = rotatedTension.xy;
  values[4] = rotatedTension.yz;
  values[5] = rotatedTension.xz;
  values[6] = values[7] = values[8] = 0;
}

template <>
void ElasticSystemCommon<Space2>::FixedBoundaryInfoFunctor::SetVelocity(
  const Vector& externalVelocity, Scalar* values)
{
  values[3] = externalVelocity.x;
  values[4] = externalVelocity.y;
}

template <>
void ElasticSystemCommon<Space3>::FixedBoundaryInfoFunctor::SetVelocity(
  const Vector& externalVelocity, Scalar* values)
{
  values[6] = externalVelocity.x;
  values[7] = externalVelocity.y;
  values[8] = externalVelocity.z;
}


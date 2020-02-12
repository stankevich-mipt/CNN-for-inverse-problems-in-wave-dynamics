#pragma once

#include "../../../../3rdparty/tinyxml/tinyxml.h"
#include "../../../../3rdparty/tinyxml/tinystr.h"
#include "ParserUtil.h"

#include <string>
#include <limits>

template<typename Space>
struct SolverSettings
{
  SPACE_TYPEDEFS

  SolverSettings(): tolerance(0.0), maxScale(1.0), maxTimeStep(1.0),
    tensionErrorMult(1.0), velocityErrorMult(1.0), positionErrorMult(1.0),
    velocityDimensionlessMult(1.0), tensionDimensionlessMult(1.0),
    allowMovement(false), allowPlasticity(false), 
    allowDiscreteDestruction(false), allowContinuousDestruction(false),
    integrator("Euler"), precision("double"), polynomialsOrder(1), hierarchyLevelsCount(1), damping(0.0),
    updateCollisionInfoPeriod(1)
  {
  }

  Scalar tolerance;
  Scalar maxScale;
  Scalar maxTimeStep;

  Scalar tensionErrorMult;
  Scalar velocityErrorMult;
  Scalar positionErrorMult;
  
  Scalar velocityDimensionlessMult;
  Scalar tensionDimensionlessMult;

  bool allowMovement;
  bool allowPlasticity;
  bool allowDiscreteDestruction;
  bool allowContinuousDestruction;

  std::string integrator;
  std::string precision;
  IndexType polynomialsOrder;
  IndexType hierarchyLevelsCount;
  Scalar damping;

  IndexType updateCollisionInfoPeriod;

  struct Erosion
  {
    Erosion() : cellAspectRatio(std::numeric_limits<Scalar>::max() / 2),
      minHeightRatio(0),
      maxPlasticDeformation(1.0),
      rhoReduction(std::numeric_limits<Scalar>::max() / 2)
    {
      boxOfInterest.Set(-std::numeric_limits<Scalar>::max() / 2 * Vector::one(), std::numeric_limits<Scalar>::max() / 2 * Vector::one());
    }
    AABB   boxOfInterest;
    Scalar cellAspectRatio;
    Scalar minHeightRatio;
    Scalar maxPlasticDeformation;
    Scalar rhoReduction;
  } erosion;

  struct DynamicContactBox
  {
    DynamicContactBox()
    {
      boxPoint1 = -Scalar(0.5) * std::numeric_limits<Scalar>::max() * Vector::one();
      boxPoint2 =  Scalar(0.5) * std::numeric_limits<Scalar>::max() * Vector::one();
    }
    Vector boxPoint1;
    Vector boxPoint2;
  } dynamicContactBox;

  void Parse(TiXmlElement* solverElement);
};

template<typename Space>
void SolverSettings<Space>::Parse(TiXmlElement *solverElement)
{
  if (ParseString(solverElement, "integrator", &integrator) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"integrator\" attribute in Solver section\n";
  }
  
  if (ParseUnsigned(solverElement, "hierarchyLevelsCount", &hierarchyLevelsCount) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"hierarchyLevelsCount\" attribute in Solver section\n";
  }

  if (ParseString(solverElement, "precision", &precision) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"precision\" attribute in Solver section\n";
  }

  if (ParseUnsigned(solverElement, "polynomialsOrder", &polynomialsOrder) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"polynomialsOrder\" attribute in Solver section\n";
  }

  if (ParseScalar(solverElement, "tolerance", &tolerance) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"tolerance\" attribute in Solver section\n";
  }

  if (ParseScalar(solverElement, "maxScale", &maxScale) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"maxScale\" attribute in Solver section\n";
  }

  if (ParseScalar(solverElement, "maxTimeStep", &maxTimeStep) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"maxTimeStep\" attribute in Solver section\n";
  }

  if (ParseBool(solverElement, "allowMovement", &allowMovement) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"allowMovement\" attribute in Solver section\n";
  }

  if (ParseBool(solverElement, "allowPlasticity", &allowPlasticity) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"allowPlasticity\" attribute in Solver section\n";
  }
      
  if (ParseBool(solverElement, "allowContinuousDestruction", &allowContinuousDestruction) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"allowContinuousDestruction\" attribute in Solver section\n";
  }

  if (ParseBool(solverElement, "allowDiscreteDestruction", &allowDiscreteDestruction) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"allowDiscreteDestruction\" attribute in Solver section\n";
  }

  if (ParseScalar(solverElement, "tensionErrorMult", &tensionErrorMult) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"tensionErrorMult\" attribute in Solver section\n";
  }

  if (ParseScalar(solverElement, "velocityErrorMult", &velocityErrorMult) != TIXML_SUCCESS)
  {
    std::cerr << "There is no required \"velocityErrorMult\" attribute in Solver section\n";
  }

  if (ParseScalar(solverElement, "positionErrorMult", &positionErrorMult) != TIXML_SUCCESS && allowMovement)
  {
    std::cerr << "\"allowMovement\" is true, however \"positionErrorMult\" is not specified\n";
  }

  if (ParseUnsigned(solverElement, "updateCollisionInfoPeriod", &updateCollisionInfoPeriod) != TIXML_SUCCESS)
  {
    std::cerr << "There is no \"updateCollisionInfoPeriod\" attribute in Solver section, using default value = 1\n";
  }

  ParseScalar(solverElement, "velocityDimensionlessMult", &velocityDimensionlessMult);  
  ParseScalar(solverElement, "tensionDimensionlessMult", &tensionDimensionlessMult);
  ParseScalar(solverElement, "damping", &damping);

  TiXmlElement* erosionElement = solverElement->FirstChildElement("Erosion");
  if (erosionElement)
  {
    ParseScalar(erosionElement, "cellAspectRatio", &erosion.cellAspectRatio);
    ParseScalar(erosionElement, "minHeightRatio",  &erosion.minHeightRatio);
    ParseScalar(erosionElement, "maxPlasticDeformation", &erosion.maxPlasticDeformation);
    ParseScalar(erosionElement, "rhoReduction", &erosion.rhoReduction);

    TiXmlElement* boxOfInterestElement = erosionElement->FirstChildElement("BoxOfInterest");
    if (boxOfInterestElement)
    {
      ParseVector(boxOfInterestElement, "boxPoint1", &erosion.boxOfInterest.boxPoint1);
      ParseVector(boxOfInterestElement, "boxPoint2", &erosion.boxOfInterest.boxPoint2);
    }
    TiXmlElement* frameOfInterestElement = erosionElement->FirstChildElement("FrameOfInterest");
    if (frameOfInterestElement)
    {
      Vector center = Vector::zeroVector();
      Vector size = Vector::one();
      Scalar margin = 0;

      ParseVector(frameOfInterestElement, "center", &center);
      ParseVector(frameOfInterestElement, "size",   &size);
      ParseScalar(frameOfInterestElement, "margin", &margin);

      erosion.boxOfInterest.boxPoint1 = center - size / Scalar(2.0) - Vector::one() * margin;
      erosion.boxOfInterest.boxPoint2 = center + size / Scalar(2.0) + Vector::one() * margin;
    }
  }

  TiXmlElement* dynamicContactBoxElement = solverElement->FirstChildElement("DynamicContactBox");
  if (dynamicContactBoxElement)
  {
    ParseVector(dynamicContactBoxElement, "boxPoint1", &dynamicContactBox.boxPoint1);
    ParseVector(dynamicContactBoxElement, "boxPoint2", &dynamicContactBox.boxPoint2);
  }

  TiXmlElement* dynamicContactFrameElement = solverElement->FirstChildElement("DynamicContactFrame");
  if (dynamicContactFrameElement)
  {
    Vector center = Vector::zeroVector();
    Vector size = Vector::one();
    Scalar margin = 0;

    ParseVector(dynamicContactFrameElement, "center", &center);
    ParseVector(dynamicContactFrameElement, "size", &size);
    ParseScalar(dynamicContactFrameElement, "margin", &margin);

    dynamicContactBox.boxPoint1 = center - size / Scalar(2.0) - Vector::one() * margin;
    dynamicContactBox.boxPoint2 = center + size / Scalar(2.0) + Vector::one() * margin;
  }
}

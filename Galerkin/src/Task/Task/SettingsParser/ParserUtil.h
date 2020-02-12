#pragma once

#include "../../../Utils/Utils.h"

template<typename Vector>
int ParseVector(TiXmlElement *element, const std::string name, Vector *vector)
{
  std::string texValue;
  element->QueryStringAttribute(name.c_str(), &texValue);
  return vector->LoadFromString(texValue) ? TIXML_SUCCESS : TIXML_WRONG_TYPE;
}

int ParseAngle(TiXmlElement *element, const std::string name, Space2::Angle* angle)
{
  std::string texValue;
  element->QueryStringAttribute(name.c_str(), &texValue);
  return sscanf(texValue.c_str(), "%lf", angle) == 1 ? TIXML_SUCCESS : TIXML_WRONG_TYPE;
}

int ParseAngle(TiXmlElement *element, const std::string name, Space3::Angle* angle)
{
  std::string texValue;
  element->QueryStringAttribute(name.c_str(), &texValue);
  return angle->LoadFromString(texValue) ? TIXML_SUCCESS : TIXML_WRONG_TYPE;
}

int ParseScalar(TiXmlElement *element, const std::string name, float *value)
{
  return element->QueryFloatAttribute(name.c_str(), value);
}

int ParseScalar(TiXmlElement *element, const std::string name, double *value)
{
  return element->QueryDoubleAttribute(name.c_str(), value);
}

int ParseBool(TiXmlElement *element, const std::string name, bool *value)
{
  return element->QueryBoolAttribute(name.c_str(), value);
}

int ParseString(TiXmlElement *element, const std::string name, std::string *value)
{
  return element->QueryStringAttribute(name.c_str(), value);
}

template<typename T>
int ParseUnsigned(TiXmlElement *element, const std::string name, T *value)
{
  int signedValue = *value;
  int res = element->QueryIntAttribute(name.c_str(), &signedValue);
  *value = signedValue;
  return res;
}

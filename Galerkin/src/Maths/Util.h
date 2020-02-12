#pragma once

#include <limits>
#include <math.h>

const double pi = 3.1415926535897932384626433;

template<typename T>
inline const T Sqr(const T &val)
{
	return val * val;
}

template<typename T>
const T& Max(const T &a, const T &b)
{
  return a > b ? a : b;
}

template<typename T>
const T& Min(const T &a, const T &b)
{
  return a < b ? a : b;
}

template<typename T>
const T Sgn(const T &c)
{
  return (c > 0.0f) ? T(1.0) : T(-1.0);
}

template <typename T>
const T Log2(const T& x)
{
  return log(x) / log(2.0);
}

template <typename T>
inline const T Abs(const T& val)
{
  return abs(val);
}

template <>
inline const double Abs(const double& val)
{
  return fabs(val);
}

template <>
inline const float Abs(const float& val)
{
  return fabs(val);
}

template <>
inline const long double Abs(const long double& val)
{
  return fabs(val);
}

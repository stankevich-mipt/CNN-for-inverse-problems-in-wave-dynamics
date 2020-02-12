#pragma once

#include "Vector2.h"
#include "Matrix2x2.h"

template<typename T>
class Tensor2
{
public:
  Tensor2()
  {
    xx = yy = xy = 0;
  }
  Tensor2(T xx, T xy, T yy)
  {
    this->xx = xx;
    this->xy = xy;
    this->yy = yy;
  }
  Tensor2(T diag)
  {
    xx = yy = diag;
    xy = 0;
  }
  Vector2<T> operator ()(const Vector2<T> &v) const
  {
    return Vector2<T>(
      xx * v.x + xy * v.y,
      xy * v.x + yy * v.y);
  }
  T GetTrace()
  {
    return xx + yy;
  }

  Tensor2& operator +=(const T &v)
  {
    this->xx += v;
    this->yy += v;
    return *this;
  }

  Tensor2& operator +=(const Tensor2& t)
  {
    this->xx += t.xx;
    this->xy += t.xy;
    this->yy += t.yy;
    return *this;
  }

  Tensor2& operator *=(const T &v)
  {
    this->xx *= v;
    this->xy *= v;
    this->yy *= v;
    return *this;
  }

  Tensor2 RotateAxes(const T& angle)
  {
    Vector2<T> n(cos(angle), sin(angle));
    return RotateAxes(n);
  }

  Tensor2 RotateAxes(const Vector2<T>& axis)
  {
    Matrix2x2<T> res(0, 0, 
                     0, 0);
    Matrix2x2<T> p(axis.x, axis.y, 
                  -axis.y, axis.x);
    Matrix2x2<T> s(xx, xy, 
                   xy, yy);

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          for (int l = 0; l < 2; ++l)
            res.data[i][j] += s.data[k][l] * p.data[i][k] * p.data[j][l];

    // assert(res.data[0][1] == res.data[1][0]);
    return Tensor2(res.data[0][0], res.data[0][1], 
                                   res.data[1][1]);
  }

  void GetEigenValues(T* eigenValues) const
  {
    T maxTangentStress = T(0.5) * (sqrt(Sqr(xx - yy) + 4 * Sqr(xy)));

    eigenValues[0] = T(0.5) * (xx + yy) + maxTangentStress;
    eigenValues[1] = T(0.5) * (xx + yy) - maxTangentStress;
  }

  void GetEigenValues(T* eigenValues, Vector2<T>* eigenVectors) const
  {
    GetEigenValues(eigenValues);
    for (int i = 0; i < 2; ++i)
      eigenVectors[i] = Vector2<T>(eigenValues[i] - yy, xy).GetNorm();

    for (int i = 0; i < 2; ++i)
    {
      if (eigenVectors[i].SquareLen() < T(1e-7))
      {
        eigenVectors[i] = eigenVectors[1 - i].GetPerpendicular();
      }
    }
  }

  T& operator[](const int i)
  {
    return *(&(xx) + i);
  }

  const T& operator[](const int i) const
  {
    return *(&(xx)+i);
  }

  T xx, xy, yy;
};

template <typename T>
Tensor2<T> TensorProduct(const Vector2<T>& v)
{
  return Tensor2<T>(
    v.x * v.x, //xx
    v.x * v.y, //xy
    v.y * v.y  //yy
  );
}

template <typename T>
Tensor2<T> SymmetricTensorProduct(const Vector2<T> &v0, const Vector2<T> &v1) {
  return Tensor2<T>(
    T(2.0) * v0.x * v1.x, //xx
    v0.x * v1.y + v0.y * v1.x, //xy
    T(2.0) * v0.y * v1.y //yy
  );
}


template<typename T>
inline const Tensor2<T> operator +(const Tensor2<T> &t0, const Tensor2<T> &t1)
{
  return Tensor2<T>(
    t0.xx + t1.xx,
    t0.xy + t1.xy,
    t0.yy + t1.yy);
}

template<typename T>
inline const Tensor2<T> operator +(const Tensor2<T> &t, const T v)
{
  Tensor2<T> res = t;
  res += v;
  return res;
}

template<typename T>
inline const Tensor2<T> operator -(const Tensor2<T> &t0, const Tensor2<T> &t1)
{
  return Tensor2<T>(
    t0.xx - t1.xx,
    t0.xy - t1.xy,
    t0.yy - t1.yy);
}

template<typename T>
inline const Tensor2<T> operator *(const Tensor2<T> &t, const T &s)
{
  return Tensor2<T>(
    t.xx * s,
    t.xy * s,
    t.yy * s);
}

template<typename T>
inline const Tensor2<T> operator *(const Tensor2<T> &t0, const Tensor2<T> &t1)
{
  return Tensor2<T>(
    t0.xx * t1.xx + t0.xy * t1.xy,
    t0.xx * t1.xy + t0.xy * t1.yy,
    t0.xy * t1.xy + t0.yy * t1.yy);
}

template <class T>
inline T DoubleConvolution(const Tensor2<T> &t0, const Tensor2<T> &t1)
{
  return 
    t0.xx * t1.xx + 
    t0.yy * t1.yy + 
    T(2.0) * t0.xy * t1.xy;
}

typedef Tensor2<float>  Tensor2f;
typedef Tensor2<double> Tensor2d;

#pragma once
#include "Vector3.h"
#include "Matrix3x3.h"
#include <assert.h>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

template<typename T>
class Tensor3
{
public:
  Tensor3()
  {
    xx = xy = xz = yy = yz = zz = 0;
  }
  Tensor3(T xx, T xy, T xz, T yy, T yz, T zz)
  {
    this->xx = xx;
    this->xy = xy;
    this->xz = xz;
    this->yy = yy;
    this->yz = yz;
    this->zz = zz;
  }
  Tensor3(T diag)
  {
    xx = yy = zz = diag;
    xy = xz = yz = 0;
  }
  Vector3<T> operator ()(const Vector3<T> &v) const
  {
    return Vector3<T>(
      xx * v.x + xy * v.y + xz * v.z,
      xy * v.x + yy * v.y + yz * v.z,
      xz * v.x + yz * v.y + zz * v.z);
  }
    
  T GetTrace() const
  {
    return xx + yy + zz;
  }

  T GetSecondInv() const
  {
    return xx * yy + xx * zz + yy * zz - Sqr(xy) - Sqr(xz) - Sqr(yz);
  }

  T GetThirdInv() const
  {
    return xx * yy * zz + 2 * xy * xz * yz 
      - xx * Sqr(yz) - yy * Sqr(xz) - zz * Sqr(xy);
  }

  T& operator[](const int i)
  {
    return *(&(xx) + i);
  }

  const T& operator[](const int i) const
  {
    return *(&(xx)+i);
  }

  Tensor3& operator +=(const T &v)
  {
    xx += v;
    yy += v;
    zz += v;
    return *this;
  }

  Tensor3& operator +=(const Tensor3& t)
  {
    xx += t.xx;
    xy += t.xy;
    xz += t.xz;
    yy += t.yy;
    yz += t.yz;
    zz += t.zz;
    return *this;
  }

  Tensor3& operator *=(const T &v)
  {
    this->xx *= v;
    this->xy *= v;
    this->xz *= v;
    this->yy *= v;
    this->yz *= v;
    this->zz *= v;
    return *this;
  }

  Tensor3 RotateAxes(const Vector3<T>& newX, const Vector3<T>& newY, const Vector3<T>& newZ)
  {
    Matrix3x3<T> res(0, 0, 0, 
                     0, 0, 0, 
                     0, 0, 0);
    Matrix3x3<T> p(newX.x, newX.y, newX.z, 
                   newY.x, newY.y, newY.z, 
                   newZ.x, newZ.y, newZ.z);
    Matrix3x3<T> s(xx, xy, xz, 
                   xy, yy, yz, 
                   xz, yz, zz);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l)
            res.data[i][j] += s.data[k][l] * p.data[i][k] * p.data[j][l];

    /*
    assert(res.data[0][1] == res.data[1][0] && 
           res.data[0][2] == res.data[2][0] && 
           res.data[1][2] == res.data[2][1]); */

    return Tensor3(res.data[0][0], res.data[0][1], res.data[0][2],
                                   res.data[1][1], res.data[1][2],
                                                   res.data[2][2]);
  }

  void GetEigenValues(T* eigenValues) const
  {
    // from http://www.continuummechanics.org/cm/principalstress.html
    T i1 = GetTrace();
    T i2 = GetSecondInv();
    T i3 = GetThirdInv();

    T q = (3 * i2 - Sqr(i1)) / 9;
    T r = (2 * pow(i1, 3) - 9 * i1 * i2 + 27 * i3) / 54;

    T theta = acos(r / pow(fabs(q), 1.5) );


    Vector3<T> res = T(2.0) * sqrt( fabs(q) ) * Vector3<T>(cos(theta / 3), cos((theta + 2 * pi) / 3), cos((theta + 4 * pi) / 3)) +
    i1 / T(3.0) * Vector3<T>::one();

    eigenValues[0] = res.x;
    eigenValues[1] = res.y;
    eigenValues[2] = res.z;
  }

  void GetEigenValues(T* eigenValues, Vector3<T>* eigenVectors) const
  {
    typedef Eigen::Matrix<T, 3, 3> Matrix;
    Matrix sigma;

    sigma << 
      xx, xy, xz,
      xy, yy, yz,
      xz, yz, zz;
    
    Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(sigma);
    
    assert(eigensolver.info() == Eigen::Success);
    
    for (int valueIndex = 0; valueIndex < 3; ++valueIndex)
    {
      eigenValues[valueIndex] = eigensolver.eigenvalues()[valueIndex];

      for (int componentIndex = 0; componentIndex < 3; ++componentIndex)
      {
        eigenVectors[valueIndex][componentIndex] = eigensolver.eigenvectors().col(valueIndex)[componentIndex];
      }
    }
  }

  T xx, xy, xz, yy, yz, zz;
};

template <typename T>
Tensor3<T> TensorProduct(const Vector3<T> &v)
{
  return Tensor3<T>(
    v.x * v.x, //xx
    v.x * v.y, //xy
    v.x * v.z, //xz
    v.y * v.y, //yy
    v.y * v.z, //yz
    v.z * v.z  //zz
  );
}

template <typename T>
Tensor3<T> SymmetricTensorProduct(const Vector3<T> &v0, const Vector3<T> &v1) {
  return Tensor3<T>(
    T(2.0) * v0.x * v1.x, //xx
    v0.x * v1.y + v0.y * v1.x, //xy
    v0.x * v1.z + v0.z * v1.x, //xz
    T(2.0) * v0.y * v1.y, //yy
    v0.y * v1.z + v0.z * v1.y, //yz
    T(2.0) * v0.z * v1.z  //zz
  );
}


template<typename T>
inline const Tensor3<T> operator +(const Tensor3<T> &t0, const Tensor3<T> &t1)
{
  return Tensor3<T>(
    t0.xx + t1.xx,
    t0.xy + t1.xy,
    t0.xz + t1.xz,
    t0.yy + t1.yy,
    t0.yz + t1.yz,
    t0.zz + t1.zz);
}

template<typename T>
inline const Tensor3<T> operator -(const Tensor3<T> &t0, const Tensor3<T> &t1)
{
  return Tensor3<T>(
    t0.xx - t1.xx,
    t0.xy - t1.xy,
    t0.xz - t1.xz,
    t0.yy - t1.yy,
    t0.yz - t1.yz,
    t0.zz - t1.zz);
}

template<typename T>
inline const Tensor3<T> operator *(const Tensor3<T> &t, const T &s)
{
  return Tensor3<T>(
    t.xx * s,
    t.xy * s,
    t.xz * s,
    t.yy * s,
    t.yz * s,
    t.zz * s);
}

template <class T>
inline T DoubleConvolution(const Tensor3<T> &v0, const Tensor3<T> &v1)
{
  return 
    v0.xx * v1.xx + 
    v0.yy * v1.yy + 
    v0.zz * v1.zz + 
    T(2.0) * v0.xy * v1.xy +
    T(2.0) * v0.xz * v1.xz +
    T(2.0) * v0.yz * v1.yz;
}

typedef Tensor3<float>  Tensor3f;
typedef Tensor3<double> Tensor3d;

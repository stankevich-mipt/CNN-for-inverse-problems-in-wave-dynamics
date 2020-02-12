#pragma once

#include "Util.h"
#include "Spaces.h"
#include <string>
#include <sstream>

template<typename T>
struct Vector2
{
  T x, y;

  T& operator[](const int i)
  {
    return *(&x + i);
  }

  T Get(const int i)
  {
    return i < 2 ? *(&x + i) : T(0.0);
  }

  const T Get(const int i) const
  {
    return i < 2 ? *(&x + i) : T(0.0);
  }

  inline Vector2<T>() {}
  inline Vector2<T>(const T& x, const T& y): x(x), y(y) { }

  template <typename U>
  Vector2<T>(const Vector2<U>& other): x(T(other.x)), y(T(other.y)) { }

  inline Vector2(const T& value): x(value), y(value)
  {}

  inline bool LoadFromString(const std::string& s)
  {
    std::stringstream stream(s);
    stream >> x >> y;
    return stream.good();
  }

  inline T Len() const
  {
    return sqrt(x * x + y * y);
  }
  
  inline T SquareLen() const
  {
    return x * x + y * y;
  }
  
  inline void Normalize() //Normalize itself
  {
    T l = sqrt(x * x + y * y);
    if (fabs(l) > T(1e-8))
    {
      T k = T(1.0) / l;
      x *= k;
      y *= k;
    }
  }
  
  inline void Invert()
  {
    x = -x;
    y = -y;
  }
  
  void Rotate(const T angle)
  {
    Vector2 self = *this;
    Vector2 x = self;
    Vector2 y = Vector2(-x.y, x.x);
    Vector2 delta = x * cos(angle) + y * sin(angle) - x;
    self += delta;
    *this = self;
  }
  inline Vector2<T> GetNorm() const
  {
    T l = sqrt(x * x + y * y);
    if (fabs(l) > T(1e-8))
    {
      T k = T(1.0) / l;
      return Vector2<T>(x * k, y * k);
    }else
    {
      return Vector2<T>(0, 0);
    }
  }

  template <typename U>
  inline Vector2<T> ComponentMultiply(const Vector2<U>& v) const
  {
    return Vector2<T>(x * v.x, y * v.y);
  }

  template <typename U>
  inline Vector2<T> ComponentDivide(const Vector2<U>& v) const
  {
    return Vector2<T>(x / v.x, y / v.y);
  }

  inline Vector2<T> ComponentAbs() const
  {
    return Vector2<T>(Abs(x), Abs(y));
  }

  inline Vector2<T> operator-() const
  {
    return Vector2<T>(-x, -y);
  }
  void Decrease(T val)
  {
    if(SquareLen() > val * val)
    {
      T len = Len();
      T scale = (len - val) / len;
      x *= scale;
      y *= scale;
    }else
    {
      x = 0.0f;
      y = 0.0f;
    }
  }
  inline Vector2<T> &operator *=(const T &i)			//mul number
  {
    x *= i;
    y *= i;
    return *this;
  }
  inline Vector2<T> &operator /=(const T &i)			//div number
  {
    T inv = T(1.0) / i;
    x *= inv;
    y *= inv;
    return *this;
  }
  inline Vector2<T> &operator +=(const Vector2<T> &V)		//plus Vector3d
  {
    x += V.x;
    y += V.y;
    return *this;
  }
  inline Vector2<T> &operator-=(const Vector2<T> &V)		//plus Vector3d
  {
    x -= V.x;
    y -= V.y;
    return *this;
  }
  inline Vector2<T> &operator--()
  {
    x = -x;
    y = -y;
    return *this;
  }

  inline Vector2<T> operator+(const Vector2<T> &V) const
  {
    //return *((Vector3d*)&_mm_add_ps(_mm_load_ps(&x), _mm_load_ps(&V.x)));
    return Vector2<T>(x + V.x, y + V.y);
  }
  inline Vector2<T> operator-(const Vector2<T> &V) const
  {
    return Vector2<T>(x - V.x, y - V.y);
  }
  inline T operator*(const Vector2<T> &V) const
  {
    return x * V.x + y * V.y;
  }
  inline Vector2<T> operator*(const T &d) const
  {
    return Vector2<T>(x * d, y * d);
  }
  Vector2<T> GetPerpendicular() const
  {
    return Vector2<T>(-y, x);
  }

  template<typename SomeType>
  inline Vector2<T>& operator=(const Vector2<SomeType>& v)
  {
    this->x = v.x;
    this->y = v.y;
    return *this;
  }

  static const Vector2<T> zeroVector()
  {
    return Vector2<T>(0, 0);
  }

  static const Vector2<T> zero()
  {
    return zeroVector();
  }

  static const Vector2<T> one()
  {
    return Vector2<T>(T(1.0), T(1.0));
  }

  static const Vector2<T> xAxis()
  {
    return Vector2<T>(T(1.0), 0);
  }

  static const Vector2<T> yAxis()
  {
    return Vector2<T>(0, T(1.0));
  }

  std::string ToString() const
  {
    std::stringstream s;
    s << x << " " << y;
    return s.str();
  }

  T GetVolume() const
  {
    return x * y;
  }
};

template<typename T>
inline Vector2<T> operator*(const T &d, const Vector2<T> &V)
{
  return Vector2<T>(V.x * d, V.y * d);
}

template<typename T>
inline Vector2<T> operator/(const Vector2<T> &V, const T &d)
{
  T invd;
  if(fabs(d) > T(1e-8)) invd = T(1.0) / d; else invd = std::numeric_limits<T>::infinity();
  return Vector2<T>(V.x * invd, V.y * invd);
}

template<typename T>
inline T operator^(const Vector2<T> &v1, const Vector2<T> &v2)
{
    return v1.x * v2.y - v1.y * v2.x;
}

typedef Vector2<float>  Vector2f;
typedef Vector2<double> Vector2d;


const static   Vector2f zeroVector2f = Vector2f(0, 0);
const static   Vector2f xAxis2f      = Vector2f(1, 0);
const static   Vector2f yAxis2f      = Vector2f(0, 1);

const static   Vector2d zeroVector2d = Vector2d(0, 0);
const static   Vector2d xAxis2d      = Vector2d(1, 0);
const static   Vector2d yAxis2d      = Vector2d(0, 1);

template<typename T>
bool TwoLinesCrossing(const Vector2<T> &p1,const Vector2<T> &p2, const Vector2<T> &t1, const Vector2<T> &t2, Vector2<T> &p0)
{
  Vector2<T> v1, v2;
  T k1, k2;
  v1 = p2 - p1;
  v2 = t2 - t1;
  T invmul;
  T mul = v1 ^ v2;
  if(fabs(mul) > T(1e-5))
  {
    invmul = 1.0f / (v1 ^ v2);
    k2 = ((t1 ^ v1) - (p1 ^ v1)) * invmul;/*p1.x * v1.y - p1.y * v1.x + t1.y * v1.x - t1.x * v1.y*/
    k1 = ((t1 ^ v2) - (p1 ^ v2)) * invmul;
    p0 = p1 + v1 * k1;
//		Vector p02 = t1 + v2 * k2;
    return((k1 > 0.0f) && (k1 < 1.0f) && (k2 > 0.0f) && (k2 < 1.0f));
  }else
  {
    return 0;
  }
  p0 = t1 + (t1 - t2); //100% bad point
  return 0;
}


template<typename T>
bool ProjectPointAgainstLine(const Vector2<T> &t1, const Vector2<T> &t2, const Vector2<T> &p, Vector2<T> &p0, T &signOfSide)
{
  Vector2<T> v1 = p - t1;
  Vector2<T> v2 = t2 - t1;
  signOfSide = Sgn(v1 ^ v2);
  p0 = t1 + v2 * ((v1 * v2) / v2.SquareLen());
  if((v1 * v2 >= 0.0f) && ((v1 * v2) / (v2.SquareLen()) <= 1.0f))
  {
    return 1;
  }
    else
  {
    return 0;
  }
}

template<typename T>
bool ProjectPointAgainstLine(const Vector2<T> &t1, const Vector2<T> &t2, const Vector2<T> &p, Vector2<T> &p0)
{
  T signOfSide;
  return ProjectPointAgainstLine(t1, t2, p, p0, signOfSide);
}

template <typename T>
T PointToSegmentDistanse(const Vector2<T> &t1, const Vector2<T> &t2, const Vector2<T> &p)
{
  Vector2<T> p0;
  T signOfSide;
  if (ProjectPointAgainstLine(t1, t2, p, p0, signOfSide)) 
  {
    return Vector2<T>(p.x - p0.x, p.y - p0.y).Len();
  } else 
  {
    return Min(Vector2<T>(p.x - t1.x, p.y - t1.y).Len(),
               Vector2<T>(p.x - t2.x, p.y - t2.y).Len());
  }
}

template<typename T>
void ProjectPointAgainstLine(const Vector2<T> &point, const Vector2<T> &planePoint, const Vector2<T> &planeNormal, const Vector2<T> &projectionDirection,
  Vector2<T> &projectedPoint)
{
  T mult = T(1.0) / (projectionDirection * planeNormal);
  projectedPoint = point + projectionDirection * ((planePoint * planeNormal) - (point * planeNormal)) * mult;
}

template<typename T>
void ProjectPointAgainstLine_noreturn(const Vector2<T> &t1, const Vector2<T> &t2, const Vector2<T> &p, Vector2<T> &p0, T &signOfSide)
{
  Vector2<T> v1 = p - t1;
  Vector2<T> v2 = t2 - t1;
  signOfSide = Sgn(v1 ^ v2);
  p0 = t1 + v2 * ((v1 * v2) / v2.SquareLen());
}

template<typename T>
bool PointInCellEx(const Vector2<T> points[3], Vector2<T> testPoint, T eps = 0)
{
  typedef T Scalar;
  Scalar side0 = ((points[1] - points[0]) ^ (testPoint - points[0]));
  Scalar side1 = ((points[2] - points[1]) ^ (testPoint - points[1]));
  Scalar side2 = ((points[0] - points[2]) ^ (testPoint - points[2]));

  if (side0 >= -eps && side1 >= -eps && side2 >= -eps) return true;
  if (side0 <=  eps && side1 <=  eps && side2 <=  eps) return true;
  return false;
}

template <typename T>
bool PointInCell(const Vector2<T> points[3], const Vector2<T>& testPoint)
{
  typedef T Scalar;

  //Scalar //eps = 0;//-1e-4;//std::numeric_limits<float>::epsilon();//Scalar(1e-9);
  Scalar eps = std::numeric_limits<Scalar>::epsilon();

  Scalar side0 = ((points[1] - points[0]) ^ (testPoint - points[0]));
  Scalar side1 = ((points[2] - points[1]) ^ (testPoint - points[1]));
  Scalar side2 = ((points[0] - points[2]) ^ (testPoint - points[2]));

  if (side0 >= -eps && side1 >= -eps && side2 >= -eps) return true;
  if (side0 <=  eps && side1 <=  eps && side2 <=  eps) return true;
  return false;

  /*Scalar eps = std::numeric_limits<float>::epsilon();//Scalar(1e-9);
  Scalar side012 =  mixed_product(points[1] - points[0], points[2] - points[0], testPoint - points[0]) *
                    mixed_product(points[1] - points[0], points[2] - points[0], points[3] - points[0]);
  if(side012 < -eps) return 0;

  Scalar side123 =  mixed_product(points[1] - points[2], points[3] - points[2], testPoint - points[2]) *
                    mixed_product(points[1] - points[2], points[3] - points[2], points[0] - points[2]);
  if(side123 < -eps) return 0;

  Scalar side230 =  mixed_product(points[2] - points[3], points[0] - points[3], testPoint - points[3]) *
                    mixed_product(points[2] - points[3], points[0] - points[3], points[1] - points[3]);
  if(side230 < -eps) return 0;

  Scalar side013 =  mixed_product(points[0] - points[1], points[3] - points[1], testPoint - points[1]) *
                    mixed_product(points[0] - points[1], points[3] - points[1], points[2] - points[1]);
  if(side013 < -eps) return 0;

  return 1;*/


  /*Scalar side1 = mixed_product(points[2] - points[0], points[3] - points[0], testPoint - points[0]);
  Scalar side2 = mixed_product(points[3] - points[0], points[1] - points[0], testPoint - points[0]);
  Scalar side3 = mixed_product(points[3] - points[1], points[2] - points[1], testPoint - points[1]);
  if (side0 * side1 < 0) return 0;
  if (side1 * side2 < 0) return 0;
  if (side2 * side3 < 0) return 0;

  return 1;*/
}

template<typename T>
T TriHeight(const Vector2<T>& v0, const Vector2<T>& v1)
{
  Vector2<T> t = v1 - v0;
  return (v0 - t * (v0 * t) / t.SquareLen()).Len();
}

template<typename T>
T TriHeight(Vector2<T> vectors[2])
{
  return TriHeight(vectors[0], vectors[1]);
}

template<typename Vector2, typename Scalar>
bool CellsIntersect(
  const Vector2 *points0,
  const Vector2 *points1)
{
  for(int side = 0; side < 2; side++)
  {
    const Vector2 *vertices0;
    const Vector2 *vertices1;
    switch(side)
    {
      case 0:
      {
        vertices0 = points0;
        vertices1 = points1;
      }break;
      case 1:
      {
        vertices0 = points1;
        vertices1 = points0;
      }break;
    }
    for(int edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      Vector2 edgePoints[2];
      edgePoints[0] = vertices0[edgeNumber];
      edgePoints[1] = vertices0[(edgeNumber + 1) % 3];
      Vector2 extraPoint = vertices0[(edgeNumber + 2) % 3];

      Vector2 tangent = edgePoints[1] - edgePoints[0];
      Vector2 norm = Vector2(-tangent.y, tangent.x);

      Scalar refProjection = (extraPoint - edgePoints[0]) * norm;

      Scalar projections[3];
      for(int pointNumber = 0; pointNumber < 3; pointNumber++)
      {
        projections[pointNumber] = (vertices1[pointNumber] - edgePoints[0]) * norm;
      }
      if(projections[0] > 0 && projections[1] > 0 && projections[2] > 0 && refProjection <= 0) return 0;
      if(projections[0] < 0 && projections[1] < 0 && projections[2] < 0 && refProjection >= 0) return 0;
    }
  }
  return 1;
}

#pragma once

#include "Vector3.h"

template<typename T>
class AABB3
{
private:
  bool empty;
public:
  AABB3()
  {
    empty = 1;
  }
  AABB3(const Vector3<T> &_boxPoint1, const Vector3<T> &_boxPoint2)
  {
    empty = 1;
    Set(_boxPoint1, _boxPoint2);
  }
  bool Intersects(const AABB3<T> &aabb) const
  {
    if(empty || aabb.empty) return 0;
    /*float min1 =      boxPoint1.x;
    float max1 =      boxPoint2.x;
    float min2 = aabb.boxPoint1.x;
    float max2 = aabb.boxPoint2.x;*/
    if((boxPoint1.x > aabb.boxPoint2.x) || (aabb.boxPoint1.x > boxPoint2.x)) return 0;
    if((boxPoint1.y > aabb.boxPoint2.y) || (aabb.boxPoint1.y > boxPoint2.y)) return 0;
    if((boxPoint1.z > aabb.boxPoint2.z) || (aabb.boxPoint1.z > boxPoint2.z)) return 0;
    return 1;
  }	
  bool Includes(const Vector3<T> &point) const
  {
    if(empty) return 0;
    if ((point.x < boxPoint1.x) || (point.x > boxPoint2.x) ||
      (point.y < boxPoint1.y) || (point.y > boxPoint2.y) ||
      (point.z < boxPoint1.z) || (point.z > boxPoint2.z))
    {
      return 0;
    }
    return 1;
  }
	bool Includes(const AABB3<T> &aabb) const
	{
    return Includes(aabb.boxPoint1) && Includes(aabb.boxPoint2);
	}

  void Set(const Vector3<T> &_boxPoint1, const Vector3<T> &_boxPoint2)
  {
    empty = 0;
    boxPoint1 = _boxPoint1;
    boxPoint2 = _boxPoint2;
  }
  void Reset()
  {
    empty = 1;
    boxPoint1 = Vector3<T>::zeroVector();
    boxPoint2 = Vector3<T>::zeroVector();
  }
  void Expand(const Vector3<T> &additionalPoint)
  {
    if(empty)
    {
      empty = 0;
      boxPoint1 = additionalPoint;
      boxPoint2 = additionalPoint;
    }else
    {
      boxPoint1.x = Min(boxPoint1.x, additionalPoint.x);
      boxPoint1.y = Min(boxPoint1.y, additionalPoint.y);
      boxPoint1.z = Min(boxPoint1.z, additionalPoint.z);

      boxPoint2.x = Max(boxPoint2.x, additionalPoint.x);
      boxPoint2.y = Max(boxPoint2.y, additionalPoint.y);
      boxPoint2.z = Max(boxPoint2.z, additionalPoint.z);
    }
  }
  void Expand(const AABB3<T> &internalAABB)
  {
    Expand(internalAABB.boxPoint1);
    Expand(internalAABB.boxPoint2);
  }
  void Expand(const T mult)
  {
    Vector3<T> size = boxPoint2 - boxPoint1;
    boxPoint1 -= size * (mult - T(1.0)) * T(0.5);
    boxPoint2 += size * (mult - T(1.0)) * T(0.5);
  }
  T Square() const
  {
    Vector3<T> size = boxPoint2 - boxPoint1;
    return (size.x * size.y + size.x * size.z + size.y * size.z) * T(2.0);
  }
  T Volume() const
  {
    Vector3<T> size = boxPoint2 - boxPoint1;
    return size.x * size.y * size.z;
  }
  Vector3<T> boxPoint1, boxPoint2;
};

typedef AABB3<float>	AABB3f;
typedef AABB3<double>	AABB3d;

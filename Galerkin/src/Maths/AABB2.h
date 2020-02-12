#pragma once

#include "Vector2.h"

template<typename T>
class AABB2
{
private:
  bool empty;
public:
  AABB2()
  {
    empty = true;
  }
  AABB2(const Vector2<T> &_boxPoint1, const Vector2<T> &_boxPoint2)
  {
    empty = true;
    Set(_boxPoint1, _boxPoint2);
  }

  bool Intersects(const AABB2<T> &aabb) const
  {
    if(empty || aabb.empty) return 0;
    /*float min1 =      boxPoint1.x;
    float max1 =      boxPoint2.x;
    float min2 = aabb.boxPoint1.x;
    float max2 = aabb.boxPoint2.x;*/
    if((boxPoint1.x > aabb.boxPoint2.x) || (aabb.boxPoint1.x > boxPoint2.x)) return 0;
    if((boxPoint1.y > aabb.boxPoint2.y) || (aabb.boxPoint1.y > boxPoint2.y)) return 0;
    return 1;
  }	
  bool Includes(const Vector2<T> &point) const
  {
    if(empty) return 0;
    if ((point.x < boxPoint1.x) || (point.x > boxPoint2.x) ||
      (point.y < boxPoint1.y) || (point.y > boxPoint2.y))
    {
      return 0;
    }
    return 1;
  }
  bool Includes(const AABB2<T> &aabb) const
  {
    return Includes(aabb.boxPoint1) && Includes(aabb.boxPoint2);
  }
  void Set(const Vector2<T> &_boxPoint1, const Vector2<T> &_boxPoint2)
  {
    empty = 0;
    boxPoint1 = _boxPoint1;
    boxPoint2 = _boxPoint2;
  }
  void Reset()
  {
    empty = true;
    boxPoint1 = Vector2<T>::zeroVector();
    boxPoint2 = Vector2<T>::zeroVector();
  }
  void Expand(const Vector2<T> &additionalPoint)
  {
    if(empty)
    {
      empty = false;
      boxPoint1 = additionalPoint;
      boxPoint2 = additionalPoint;
    }else
    {
      boxPoint1.x = Min(boxPoint1.x, additionalPoint.x);
      boxPoint1.y = Min(boxPoint1.y, additionalPoint.y);

      boxPoint2.x = Max(boxPoint2.x, additionalPoint.x);
      boxPoint2.y = Max(boxPoint2.y, additionalPoint.y);
    }
  }
  void Expand(const AABB2<T> &internalAABB)
  {
    Expand(internalAABB.boxPoint1);
    Expand(internalAABB.boxPoint2);
  }
  void Expand(const T mult)
  {
    Vector2<T> size = boxPoint2 - boxPoint1;
    boxPoint1 -= size * (mult - T(1.0)) * T(0.5);
    boxPoint2 += size * (mult - T(1.0)) * T(0.5);
  }
  T Square() const //perimeter actually, used for AABBTree balancing
  {
    Vector2<T> size = boxPoint2 - boxPoint1;
    return (size.x + size.y) * T(2.0);
  }
  T Volume() const
  {
    Vector2<T> size = boxPoint2 - boxPoint1;
    return size.x * size.y;
  }
  Vector2<T> boxPoint1, boxPoint2;
};

typedef AABB2<float>  AABB2f;
typedef AABB2<double> AABB2d;

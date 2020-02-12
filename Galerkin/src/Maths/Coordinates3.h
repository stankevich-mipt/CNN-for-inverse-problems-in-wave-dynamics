#pragma once
#include "Vector3.h"

template<typename T>
class Coordinates3
{
public:
  Vector3<T> xVector, yVector, zVector;
  Vector3<T> pos;
  Coordinates3(){}
  Coordinates3(const Vector3<T> &_pos, const Vector3<T> &_xVector, const Vector3<T> &_yVector, const Vector3<T> &_zVector)
  {
    xVector = _xVector;
    yVector = _yVector;
    zVector = _zVector;

    pos = _pos;
  }

  /*static const Coordinates Identity()
  {
    return Coordinates(zeroVector3<T>, xAxis, yAxis, zAxis);
  }*/

  const Vector3<T> GetPointRelativePos(const Vector3<T> &globalPoint) const
  {
    Vector3<T> delta = globalPoint - pos;
    return Vector3<T>(delta * xVector, 
                      delta * yVector, 
                      delta * zVector);
  }

  const Vector3<T> GetAxisRelativeOrientation(const Vector3<T> &globalAxis) const 
  {
    Vector3<T> delta = globalAxis;
    return Vector3<T>(delta * xVector, 
                      delta * yVector, 
                      delta * zVector);
  }
  const Vector3<T> GetPointGlobalPos  (const Vector3<T> &relativePoint) const 
  {
    return pos + xVector * relativePoint.x + yVector * relativePoint.y + zVector * relativePoint.z;
  }
  const Vector3<T> GetAxisGlobalOrientation  (const Vector3<T> &relativeAxis) const 
  {
    return xVector * relativeAxis.x + yVector * relativeAxis.y + zVector * relativeAxis.z;
  }
  const Coordinates3<T> GetGlobalCoordinates(const Coordinates3<T> &localCoordinates)
  {
    return Coordinates3<T>( GetPointGlobalPos(localCoordinates.pos),
                            GetAxisGlobalOrientation(localCoordinates.xVector),
                            GetAxisGlobalOrientation(localCoordinates.yVector),
                            GetAxisGlobalOrientation(localCoordinates.zVector));
  }

  const Coordinates3<T> GetLocalCoordinates(const Coordinates3<T> &globalCoordinates)
  {
    return Coordinates3<T>( GetPointRelativePos(globalCoordinates.pos),
                            GetAxisRelativeOrientation(globalCoordinates.xVector),
                            GetAxisRelativeOrientation(globalCoordinates.yVector),
                            GetAxisRelativeOrientation(globalCoordinates.zVector));
  }

  void Rotate(Vector3<T> point, Vector3<T> axis)
  {
    T angle = axis.Len();
    if(angle == 0) return;
    Vector3f normAxis = axis / angle;
    Vector3f delta = point - this->pos;
    delta.Rotate(normAxis, angle);
    pos = point - delta;

    xVector.Rotate(normAxis, angle);
    yVector.Rotate(normAxis, angle);
    zVector.Rotate(normAxis, angle);
  }

  void Rotate(Vector3<T> point, Vector3<T> initialAxis, Vector3<T> dstAxis)
  {
    Vector3f delta = point - this->pos;
    delta.Rotate(initialAxis, dstAxis);
    pos = point - delta;

    xVector.Rotate(initialAxis, dstAxis);
    yVector.Rotate(initialAxis, dstAxis);
    zVector.Rotate(initialAxis, dstAxis);
  }

  void Identity()
  {
    pos = Vector3<T>(0, 0, 0);
    xVector = Vector3<T>(1, 0, 0);
    yVector = Vector3<T>(0, 1, 0);
    zVector = Vector3<T>(0, 0, 1);
  }

  void Normalize()
  {
    Vector3f oldXVector = xVector;
    Vector3f oldYVector = yVector;
    Vector3f oldZVector = zVector;
    xVector = (yVector ^ zVector).GetNorm();
    yVector = (zVector ^ xVector).GetNorm();
    zVector = (xVector ^ yVector).GetNorm();

    if(oldXVector * xVector < 0) xVector.Invert();
    if(oldYVector * yVector < 0) yVector.Invert();
    if(oldZVector * zVector < 0) zVector.Invert();
  }

  const Vector3<T>& GetGlobalxVector(){return xVector;}
  const Vector3<T>& GetGlobalyVector(){return yVector;}
  const Vector3<T>& GetGlobalzVector(){return zVector;}
  const Vector3<T>& GetGlobalPos(){return pos;}

  static const Coordinates3<T> defCoords()
  {
    return Coordinates3(Vector3<T>::zeroVector(), Vector3<T>::xAxis(), Vector3<T>::yAxis(), Vector3<T>::zAxis());
  }
};

typedef Coordinates3<float>  Coordinates3f;
typedef Coordinates3<double> Coordinates3d;

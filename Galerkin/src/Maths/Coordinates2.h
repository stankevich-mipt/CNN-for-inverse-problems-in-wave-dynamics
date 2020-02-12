#pragma once
#include "Vector2.h"

template<typename T>
class Coordinates2
{
public:
	T angle;
	Vector2<T> xVector, yVector;
	Vector2<T> pos;
	Coordinates2(){}
	Coordinates2(const Vector2<T> &_pos, const T _angle)
	{
		angle = _angle;
		xVector = Vector2<T>(cos(angle), sin(angle));
		yVector = Vector2<T>(cos(angle + T(pi) / T(2.0)), sin(angle + T(pi) / T(2.0)));

		pos = _pos;
	}

	const Vector2<T> GetPointRelativePos(const Vector2<T> &globalPoint) const
	{
		Vector2<T> delta = globalPoint - pos;
		return Vector2<T>(delta * xVector, 
						delta * yVector);
	}

	const Vector2<T> GetAxisRelativeOrientation(const Vector2<T> &globalAxis) const
	{
		Vector2<T> delta = globalAxis;
		return Vector2<T>(delta * xVector, 
						delta * yVector);
	}

	const Vector2<T> GetPointGlobalPos  (const Vector2<T> &relativePoint) const
	{
		return pos + xVector * relativePoint.x + yVector * relativePoint.y;
	}

	const Vector2<T> GetAxisGlobalOrientation  (const Vector2<T> &relativeAxis) const
	{
		return xVector * relativeAxis.x + yVector * relativeAxis.y;
	}

	const Coordinates2 GetGlobalCoordinates(const Coordinates2 &localCoordinates)
	{
		return Coordinates2(GetPointGlobalPos(localCoordinates.pos), localCoordinates.angle + angle);
	}

	const Coordinates2 GetLocalCoordinates(const Coordinates2 &globalCoordinates)
	{
		return Coordinates2(GetPointRelativePos(globalCoordinates.pos), globalCoordinates.angle - angle);
	}

	void Identity()
	{
		angle = 0.0f;
		xVector = Vector2<T>(1.0f, 0.0f);
		yVector = Vector2<T>(0.0f, 1.0f);
    pos = Vector2<T>::zeroVector();
	}

	void SetRotation(const T &_angle)
	{
		angle = _angle;
		xVector = Vector2<T>(cos(angle), sin(angle));
		yVector = Vector2<T>(cos(angle + pi / 2.0), sin(angle + pi / 2.0));
	}

  static const Coordinates2<T> defCoords()
  {
    Coordinates2<T> coords;
    coords.pos = Vector2<T>::zeroVector();
    coords.xVector = Vector2<T>::xAxis();
    coords.yVector = Vector2<T>::yAxis();
    coords.angle = 0;
    return coords;
  }

};

typedef Coordinates2<float>	Coordinates2f;
typedef Coordinates2<double>	Coordinates2d;

#pragma once
#include "Vector2.h"


template<typename T>
class Matrix2x2
{
private:
	void Swap(int i1, int j1, int i2, int j2)
	{
		T tmp = data[i1][j1];
		data[i1][j1] = data[i2][j2];
		data[i2][j2] = tmp;
	}
public:
	Matrix2x2()
	{
		data[0][0] = 0; data[0][1] = 0;
		data[1][0] = 0; data[1][1] = 0;
	}
	Matrix2x2(T a11, T a12,
			  T a21, T a22)
	{
		data[0][0] = a11; data[0][1] = a12;
		data[1][0] = a21; data[1][1] = a22;
	}
	Matrix2x2(const Vector2<T> &column1, const Vector2<T> &column2)
	{
		data[0][0] = column1.x; data[0][1] = column2.x;
		data[1][0] = column1.y; data[1][1] = column2.y;
	}
	void Transpose()
	{
		Swap(1, 2, 2, 1);
	}
	inline const Matrix2x2<T> operator*(const T s) const
	{
		return Matrix2x2(data[0][0] * s, data[0][1] * s,
						 data[1][0] * s, data[1][1] * s);
	}
	inline const Matrix2x2<T> operator+(const Matrix2x2<T> m) const
	{
		return Matrix2x2(data[0][0] + m.data[0][0], data[0][1] + m.data[0][1],
						 data[1][0] + m.data[1][0], data[1][1] + m.data[1][1]);
	}
	inline const Matrix2x2<T> operator-(const Matrix2x2<T> m) const
	{
		return Matrix2x2(data[0][0] - m.data[0][0], data[0][1] - m.data[0][1],
						 data[1][0] - m.data[1][0], data[1][1] - m.data[1][1]);
	}
	const Matrix2x2<T> GetTransposed() const
	{
		Matrix2x2<T> tmp(	data[0][0], data[1][0],
							data[0][1], data[1][1]);
		return tmp;
	}
	const T GetDeterminant() const
	{
		return data[0][0] * data[1][1] - data[0][1] * data[1][0];
	}
	const Matrix2x2<T> GetInverted() const
	{
		T invdet = T(1.0) / GetDeterminant();
		Matrix2x2<T> tmp(data[1][1], -data[1][0], -data[0][1], data[0][0]);
		return tmp * invdet;
	}


	const Matrix2x2<T> GetSquareRoot(int iterations) const
	{
		Matrix2x2<T> Y_(zeroVector2d, zeroVector2d, zeroVector2d);
		Y_ = *this;
        Matrix2x2<T> Z_(Vector2<T>::xAxis(), Vector2<T>::yAxis(), Vector2<T>::zAxis());
		Matrix2x2<T> Y(zeroVector2d, zeroVector2d, zeroVector2d), Z(zeroVector2d, zeroVector2d, zeroVector2d);
		for(int i = 0; i < iterations; i++)
		{
			Y = (Y_ + Z_.GetInverted()) * T(0.5);
			Z = (Z_ + Y_.GetInverted()) * T(0.5);
			Y_ = Y;
			Z_ = Z;
		}
		return Y;
	}

	const Matrix2x2<T> GetSquareRootNewtonian(int iterations, Matrix2x2<T> warmStart = identity()) const
	{
		Matrix2x2<T> currSolution = warmStart;
		for(int i = 0; i < iterations; i++)
		{
			currSolution = currSolution - (currSolution * currSolution - *this) * (currSolution * T(2.0)).GetInverted();
		}
		return currSolution;
	}

	void Identify()
	{
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 2; j++)
				if(i == j) data[i][j] = 1; else data[i][j] = 0;
	}
	const T GetQuadNorm() const
	{
		return	sqr(data[0][0]) + sqr(data[0][1]) + sqr(data[1][0]) + sqr(data[1][1]);
	}
	inline void operator += (const Matrix2x2<T> &m)
	{
		data[0][0] += m.data[0][0]; data[0][1] += m.data[0][1];
		data[1][0] += m.data[1][0]; data[1][1] += m.data[1][1];
	}

	void operator =(const Matrix2x2<T> &m)
	{
		for(int i = 0; i < 2; i++)
		{
			for(int j = 0; j < 2; j++)
			{
				data[i][j] = m.data[i][j];
			}
		}
	}

	static const Matrix2x2 identity()
	{
		Matrix2x2 i;
		i.Identify();
		return i;
	}

	T data[2][2];
};

template<typename T>
inline const Matrix2x2<T> operator*(const Matrix2x2<T> &m1, const Matrix2x2<T> &m2)
{
	return Matrix2x2<T>
          (m1.data[0][0] * m2.data[0][0] + m1.data[0][1] * m2.data[1][0],
					 m1.data[0][0] * m2.data[0][1] + m1.data[0][1] * m2.data[1][1],
					 m1.data[1][0] * m2.data[0][0] + m1.data[1][1] * m2.data[1][0],
					 m1.data[1][0] * m2.data[0][1] + m1.data[1][1] * m2.data[1][1]);
}

template<typename T>
inline const Vector2<T> operator*(const Vector2<T> &v, const Matrix2x2<T> &m)
{
	return Vector2<T>(	v.x * m.data[0][0] + v.y * m.data[1][0],
						v.x * m.data[0][1] + v.y * m.data[1][1]);
}

template<typename T>
inline const Vector2<T> operator*(const Matrix2x2<T> &m, const Vector2<T> &v)
{
	return Vector2<T>(	v.x * m.data[0][0] + v.y * m.data[0][1],
						v.x * m.data[1][0] + v.y * m.data[1][1]);
}

typedef Matrix2x2<float>	Matrix2x2f;
typedef Matrix2x2<double> Matrix2x2d;

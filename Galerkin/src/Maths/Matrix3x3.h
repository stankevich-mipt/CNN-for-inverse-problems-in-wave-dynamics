#pragma once
#include "Vector3.h"

template<typename T>
class Matrix3x3
{
private:
	void Swap(int i1, int j1, int i2, int j2)
	{
		T tmp = data[i1][j1];
		data[i1][j1] = data[i2][j2];
		data[i2][j2] = tmp;
	}
public:
	Matrix3x3()
	{
		data[0][0] = 0; data[0][1] = 0; data[0][2] = 0; 
		data[1][0] = 0; data[1][1] = 0; data[1][2] = 0; 
		data[2][0] = 0; data[2][1] = 0; data[2][2] = 0; 
	}
	Matrix3x3(T a11, T a12, T a13, 
			  T a21, T a22, T a23, 
			  T a31, T a32, T a33)
	{
		data[0][0] = a11; data[0][1] = a12; data[0][2] = a13; 
		data[1][0] = a21; data[1][1] = a22; data[1][2] = a23; 
		data[2][0] = a31; data[2][1] = a32; data[2][2] = a33; 
	}
	Matrix3x3(const Vector3<T> &column1, const Vector3<T> &column2, const Vector3<T> &column3)
	{
		data[0][0] = column1.x; data[0][1] = column2.x; data[0][2] = column3.x; 
		data[1][0] = column1.y; data[1][1] = column2.y; data[1][2] = column3.y; 
		data[2][0] = column1.z; data[2][1] = column2.z; data[2][2] = column3.z; 
	}
	void Transpose()
	{
		Swap(1, 2, 2, 1);
		Swap(1, 3, 3, 1);
		Swap(3, 2, 2, 3);
	}
	inline const Matrix3x3<T> operator*(const T s) const
	{
		return Matrix3x3(data[0][0] * s, data[0][1] * s, data[0][2] * s,
						 data[1][0] * s, data[1][1] * s, data[1][2] * s,
						 data[2][0] * s, data[2][1] * s, data[2][2] * s);
	}
	inline const Matrix3x3<T> operator+(const Matrix3x3<T> m) const
	{
		return Matrix3x3(data[0][0] + m.data[0][0], data[0][1] + m.data[0][1], data[0][2] + m.data[0][2], 
						 data[1][0] + m.data[1][0], data[1][1] + m.data[1][1], data[1][2] + m.data[1][2], 
						 data[2][0] + m.data[2][0], data[2][1] + m.data[2][1], data[2][2] + m.data[2][2]);
	}
	inline const Matrix3x3<T> operator-(const Matrix3x3<T> m) const
	{
		return Matrix3x3(data[0][0] - m.data[0][0], data[0][1] - m.data[0][1], data[0][2] - m.data[0][2], 
						 data[1][0] - m.data[1][0], data[1][1] - m.data[1][1], data[1][2] - m.data[1][2], 
						 data[2][0] - m.data[2][0], data[2][1] - m.data[2][1], data[2][2] - m.data[2][2]);
	}
	const Matrix3x3<T> GetTransposed() const
	{
		Matrix3x3<T> tmp(	data[0][0], data[1][0], data[2][0],
							data[0][1], data[1][1], data[2][1],
							data[0][2], data[1][2], data[2][2]);
		return tmp;
	}
	const T GetDeterminant() const
	{
		return data[0][1] * data[1][2] * data[2][0] + data[0][0] * data[1][1] * data[2][2] + data[0][2] * data[1][0] * data[2][1]
			  -data[0][0] * data[2][1] * data[1][2] - data[2][2] * data[0][1] * data[1][0] - data[2][0] * data[1][1] * data[0][2];
	}
	const Matrix3x3<T> GetInverted() const
	{
		T invdet = 1.0f / GetDeterminant();
		Matrix3x3<T>     tmp(	(data[1][1] * data[2][2] - data[1][2] * data[2][1]), -(data[0][1] * data[2][2] - data[0][2] * data[2][1]),  (data[0][1] * data[1][2] - data[0][2] * data[1][1]),
						  -		(data[1][0] * data[2][2] - data[1][2] * data[2][0]),  (data[0][0] * data[2][2] - data[0][2] * data[2][0]), -(data[0][0] * data[1][2] - data[0][2] * data[1][0]),
								(data[1][0] * data[2][1] - data[1][1] * data[2][0]), -(data[0][0] * data[2][1] - data[0][1] * data[2][0]),  (data[0][0] * data[1][1] - data[0][1] * data[1][0]));
		return tmp * invdet;
	}


	const Matrix3x3<T> GetSquareRoot(int iterations) const
	{
		Matrix3x3<T> Y_(zeroVector3d, zeroVector3d, zeroVector3d);
		Y_ = *this;
        Matrix3x3<T> Z_(Vector3<T>::xAxis(), Vector3<T>::yAxis(), Vector3<T>::zAxis());
		Matrix3x3<T> Y(zeroVector3d, zeroVector3d, zeroVector3d), Z(zeroVector3d, zeroVector3d, zeroVector3d);
		for(int i = 0; i < iterations; i++)
		{
			Y = (Y_ + Z_.GetInverted()) * 0.5f;
			Z = (Z_ + Y_.GetInverted()) * 0.5f;
			Y_ = Y;
			Z_ = Z;
		}
		return Y;
	}

	const Matrix3x3<T> GetSquareRootNewtonian(int iterations, Matrix3x3<T> warmStart = identity()) const
	{
		Matrix3x3<T> currSolution = warmStart;
		for(int i = 0; i < iterations; i++)
		{
			currSolution = currSolution - (currSolution * currSolution - *this) * (currSolution * 2.0).GetInverted();
		}
		return currSolution;
	}

	void Identify()
	{
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				if(i == j) data[i][j] = 1; else data[i][j] = 0;
	}
	const T GetQuadNorm() const
	{
		return	sqr(data[0][0]) + sqr(data[0][1]) + sqr(data[0][2]) +
				sqr(data[1][0]) + sqr(data[1][1]) + sqr(data[1][2]) +
				sqr(data[2][0]) + sqr(data[2][1]) + sqr(data[2][2]);
	}
	inline void operator += (const Matrix3x3<T> &m)
	{
		data[0][0] += m.data[0][0]; data[0][1] += m.data[0][1]; data[0][2] += m.data[0][2];
		data[1][0] += m.data[1][0]; data[1][1] += m.data[1][1]; data[1][2] += m.data[1][2];
		data[2][0] += m.data[2][0]; data[2][1] += m.data[2][1]; data[2][2] += m.data[2][2];
	}

	void operator =(const Matrix3x3<T> &m)
	{
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				data[i][j] = m.data[i][j];
			}
		}
	}

	static const Matrix3x3 identity()
	{
		Matrix3x3 i;
		i.Identify();
		return i;
	}

	T data[3][3];
};

template<typename T>
inline const Matrix3x3<T> operator*(const Matrix3x3<T> &m1, const Matrix3x3<T> &m2)
{
	return Matrix3x3<T>(m1.data[0][0] * m2.data[0][0] + m1.data[0][1] * m2.data[1][0] + m1.data[0][2] * m2.data[2][0], 
					 m1.data[0][0] * m2.data[0][1] + m1.data[0][1] * m2.data[1][1] + m1.data[0][2] * m2.data[2][1],
					 m1.data[0][0] * m2.data[0][2] + m1.data[0][1] * m2.data[1][2] + m1.data[0][2] * m2.data[2][2],

					 m1.data[1][0] * m2.data[0][0] + m1.data[1][1] * m2.data[1][0] + m1.data[1][2] * m2.data[2][0], 
					 m1.data[1][0] * m2.data[0][1] + m1.data[1][1] * m2.data[1][1] + m1.data[1][2] * m2.data[2][1],
					 m1.data[1][0] * m2.data[0][2] + m1.data[1][1] * m2.data[1][2] + m1.data[1][2] * m2.data[2][2],

					 m1.data[2][0] * m2.data[0][0] + m1.data[2][1] * m2.data[1][0] + m1.data[2][2] * m2.data[2][0], 
					 m1.data[2][0] * m2.data[0][1] + m1.data[2][1] * m2.data[1][1] + m1.data[2][2] * m2.data[2][1],
					 m1.data[2][0] * m2.data[0][2] + m1.data[2][1] * m2.data[1][2] + m1.data[2][2] * m2.data[2][2]);
}

template<typename T>
inline const Vector3<T> operator*(const Vector3<T> &v, const Matrix3x3<T> &m)
{
	return Vector3<T>(	v.x * m.data[0][0] + v.y * m.data[1][0] + v.z * m.data[2][0],
						v.x * m.data[0][1] + v.y * m.data[1][1] + v.z * m.data[2][1],
						v.x * m.data[0][2] + v.y * m.data[1][2] + v.z * m.data[2][2]);
}

template<typename T>
inline const Vector3<T> operator*(const Matrix3x3<T> &m, const Vector3<T> &v)
{
	return Vector3<T>(	v.x * m.data[0][0] + v.y * m.data[0][1] + v.z * m.data[0][2],
						v.x * m.data[1][0] + v.y * m.data[1][1] + v.z * m.data[1][2],
						v.x * m.data[2][0] + v.y * m.data[2][1] + v.z * m.data[2][2]);
}

typedef Matrix3x3<float>	Matrix3x3f;
typedef Matrix3x3<double> Matrix3x3d;

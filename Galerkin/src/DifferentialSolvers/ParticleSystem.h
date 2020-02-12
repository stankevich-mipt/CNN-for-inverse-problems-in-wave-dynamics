#include <vector>
#include "DifferentialSystem.h"

struct Particle
{
  Particle(Vector3d _pos, Vector3d _velocity)
  {
    pos = _pos;
    velocity = _velocity;
    acceleration = zeroVector3d;
  }
	Particle(){}


  void AddAcceleration(Vector3d acc)
  {
    acceleration += acc;
  }

  void Push(Vector3d destPos)
	{
		Vector3d delta = destPos - pos;
		/*velocity += delta;
		pos += delta;*/
    acceleration += delta * 100.1f;
	}

	void PushWithFriction(Vector3d destPos, double frictionCoefficient, Vector3d obstacleVelocity = zeroVector3d)
	{
		Vector3d xVec = (destPos - pos);

		double depth = xVec.Len();

		Vector3d xProj = xVec * (xVec * (velocity - obstacleVelocity)) * (1.0f / xVec.SquareLen());
		Vector3d yProj = (velocity - obstacleVelocity) - xProj;

		yProj.Decrease(depth * frictionCoefficient);

		Vector3d decreasedVelocity = xProj + yProj + obstacleVelocity;

		Push(destPos);
    acceleration += (velocity - decreasedVelocity) * 0.001;
	}
	Vector3d pos, velocity;
	Vector3d acceleration;
};

struct Triangle
{
	Triangle(){}
	Triangle(int _p0, int _p1, int _p2)
	{
		p0 = _p0;
		p1 = _p1;
		p2 = _p2;
	}
	int p0, p1, p2;
};

struct Constraint
{
	Constraint(){}
	Constraint(Particle *_p1, Particle *_p2, double _stiffness)
	{
		p1 = _p1;
		p2 = _p2;
		valid = 1;
    stiffness = _stiffness;
		defLen = (p1->pos - p2->pos).Len();
	};
	void Satisfy()
	{
		if(!valid) return;
		Vector3d delta	= p1->pos - p2->pos;
		Vector3d norm	= delta.GetNorm();
		double   len		= norm * delta;

		/*double scaleFactor = 1.05f;
		if(len > defLen * scaleFactor)
		{
			defLen = len / scaleFactor;
		}
		if(len < defLen / scaleFactor)
		{
			defLen = len * scaleFactor;
		}*/


		Vector3d force = norm * (defLen - len) * stiffness;
    p1->AddAcceleration( force);
    p2->AddAcceleration(-force);
	}
	Particle *p1, *p2;
	double defLen;
  double stiffness;
	bool valid;
};


template<int _particlesCount, int _trianglesCount>
struct MeshlessConstraint
{
	MeshlessConstraint()
	{
		stiffness		= 0.1;
		compressibility = 1.0;
		maxDisplacement = 10.0;
		restitution		= 1.0;
    viscosity = 0.00f;
	}

	double TetrahedronVolume(Vector3d p1, Vector3d p2, Vector3d p3)
	{
		return (p2 ^ p1) * p3 / 6.0;
	}

	double TriangleSquare(Vector3d p1, Vector3d p2, Vector3d p3)
	{
		return ((p2 - p1) ^ (p3 - p1)).Len() / 2.0;
	}

	Vector3d TriangleNorm(Vector3d p1, Vector3d p2, Vector3d p3)
	{
		return ((p2 - p1) ^ (p3 - p1)).GetNorm();
	}

	void Create(Particle **_particles, Triangle *_triangles)
	{
		for(int i = 0; i < particlesCount; i++)
		{
			particles[i] = _particles[i];
		}
		for(int i = 0; i < trianglesCount; i++)
		{
			triangles[i] = _triangles[i];
		}

		for(int i = 0; i < particlesCount; i++)
		{
			destPos[i] = particles[i]->pos;
		}

//		deformation = identityMatrix;

		Matrix3x3d tmp(zeroVector3d, zeroVector3d, zeroVector3d);
		Vector3d p, q;
		Vector3d virtualMassCenter;

		virtualMassCenter = zeroVector3d;
		for(int i = 0; i < particlesCount; i++)
		{
			virtualMassCenter += destPos[i];
		}
		virtualMassCenter /= particlesCount;

		for(int i = 0; i < particlesCount; i++)
		{
			q = destPos[i] - virtualMassCenter;
			tmp.data[0][0] = q.x * q.x; tmp.data[0][1] = q.x * q.y; tmp.data[0][2] = q.x * q.z;
			tmp.data[1][0] = q.y * q.x; tmp.data[1][1] = q.y * q.y; tmp.data[1][2] = q.y * q.z;
			tmp.data[2][0] = q.z * q.x; tmp.data[2][1] = q.z * q.y; tmp.data[2][2] = q.z * q.z;
			Aqq += tmp;
		}

		Aqq = Aqq.GetInverted();

		initialVolume = 0.0;
		for(int i = 0; i < trianglesCount; i++)
		{
			initialVolume += TetrahedronVolume(particles[triangles[i].p0]->pos, particles[triangles[i].p1]->pos, particles[triangles[i].p2]->pos);
		}

		active = 1;
	}

	void Satisfy()
	{
		if(!active) return;
    Matrix3x3d E(Vector3d::xAxis(), Vector3d::yAxis(), Vector3d::zAxis());
		Matrix3x3d Apq(zeroVector3d, zeroVector3d, zeroVector3d);
		Matrix3x3d S(zeroVector3d, zeroVector3d, zeroVector3d);
		Matrix3x3d tmp(zeroVector3d, zeroVector3d, zeroVector3d);
		Vector3d realMassCenter = zeroVector3d;
		Vector3d virtualMassCenter = zeroVector3d;

		realMassCenter = zeroVector3d;
		for(int i = 0; i < particlesCount; i++)
		{
			realMassCenter += particles[i]->pos;
		}
		realMassCenter /= particlesCount;

		virtualMassCenter = zeroVector3d;
		for(int i = 0; i < particlesCount; i++)
		{
			virtualMassCenter += destPos[i];
		}
		virtualMassCenter /= particlesCount;

		Vector3d q;
		Vector3d p;

		for(int i = 0; i < particlesCount; i++)
		{
			p = particles[i]->pos - realMassCenter;
			q = /*deformation * */(destPos[i] - virtualMassCenter);
			tmp.data[0][0] = p.x * q.x; tmp.data[0][1] = p.x * q.y; tmp.data[0][2] = p.x * q.z;
			tmp.data[1][0] = p.y * q.x; tmp.data[1][1] = p.y * q.y; tmp.data[1][2] = p.y * q.z;
			tmp.data[2][0] = p.z * q.x; tmp.data[2][1] = p.z * q.y; tmp.data[2][2] = p.z * q.z;
			Apq += tmp;
		}

		Matrix3x3d R;

		Matrix3x3d A = Apq * Aqq;

		Matrix3x3d quadS = A.GetTransposed() * A;
		//S = quadS.GetSquareRoot(10);
		//double maxval = 1.0;
		double maxval = fabs(quadS.data[0][0]);
		if(maxval < fabs(quadS.data[0][1])) maxval = fabs(quadS.data[0][1]);
		if(maxval < fabs(quadS.data[0][2])) maxval = fabs(quadS.data[0][2]);
		if(maxval < fabs(quadS.data[1][0])) maxval = fabs(quadS.data[1][0]);
		if(maxval < fabs(quadS.data[1][1])) maxval = fabs(quadS.data[1][1]);
		if(maxval < fabs(quadS.data[1][2])) maxval = fabs(quadS.data[1][2]);
		if(maxval < fabs(quadS.data[2][0])) maxval = fabs(quadS.data[2][0]);
		if(maxval < fabs(quadS.data[2][1])) maxval = fabs(quadS.data[2][1]);
		if(maxval < fabs(quadS.data[2][2])) maxval = fabs(quadS.data[2][2]);

    Matrix3x3d warmStart = Matrix3x3d::identity();
		S = (quadS * (1.0 / maxval)).GetSquareRootNewtonian(5, warmStart) * (sqrt(maxval));
		warmStart = S;

		Matrix3x3d test = quadS + (S * S) * (-1.0);
		R = A * S.GetInverted();


		/*Matrix3x3 U = Apq;
		Matrix3x3 P;
		for( int i = 0; i < 20; ++i )
		{
			Matrix3x3 tr = U.GetTransposed() * U;
			P = ( E + tr ) * 0.5;
			U = U * P.GetInverted();
		}
		P = U.GetTransposed() * Apq;
		R = U;
		S = P;*/



		Vector3d tmpx = Vector3d(R.data[0][0], R.data[1][0], R.data[2][0]);
		Vector3d tmpy = Vector3d(R.data[0][1], R.data[1][1], R.data[2][1]);
		Vector3d tmpz = Vector3d(R.data[0][2], R.data[1][2], R.data[2][2]);


		/*Matrix3x3 resDeformation = S * deformation;
		if((tmpx ^ tmpy) * tmpz < 0.0f)
		{
			if(fabs(resDeformation.data[0][0]) < fabs(resDeformation.data[1][1]) && fabs(resDeformation.data[0][0]) < fabs(resDeformation.data[2][2]))
			{
				tmpx = -tmpx;
			}else
			if(fabs(resDeformation.data[1][1]) < fabs(resDeformation.data[0][0]) && fabs(resDeformation.data[1][1]) < fabs(resDeformation.data[2][2]))
			{
				tmpy = -tmpy;
			}else
			{
				tmpz = -tmpz;
			}
		}*/
/*		tmpz = tmpx ^ tmpy;
		tmpy = tmpz ^ tmpx;*/
		R.data[0][0] = tmpx.x; R.data[1][0] = tmpx.y, R.data[2][0] = tmpx.z;
		R.data[0][1] = tmpy.x; R.data[1][1] = tmpy.y, R.data[2][1] = tmpy.z;
		R.data[0][2] = tmpz.x; R.data[1][2] = tmpz.y, R.data[2][2] = tmpz.z;

    //R.Identify();


		/*double beta = 0.995f;
		Matrix3x3 resultMatrix = A * beta + R * (1.0f - beta);*/


		Matrix3x3d resultMatrix = R;//A * beta + R * (1.0f - beta);
		rotation = R;
		scale = S;



		//volume conservation
		if(1)
		{
			double volume = 0.0;
			for(int i = 0; i < trianglesCount; i++)
			{
				volume += TetrahedronVolume(particles[triangles[i].p0]->pos, particles[triangles[i].p1]->pos, particles[triangles[i].p2]->pos);
			}

			double surface = 0.0;
			for(int i = 0; i < trianglesCount; i++)
			{
				surface += TriangleSquare(particles[triangles[i].p0]->pos, particles[triangles[i].p1]->pos, particles[triangles[i].p2]->pos);
			}

			double linearDelta = compressibility * (initialVolume - volume) / surface;
			for(int i = 0; i < trianglesCount; i++)
			{
				Particle *p0 = particles[triangles[i].p0];
				Particle *p1 = particles[triangles[i].p1];
				Particle *p2 = particles[triangles[i].p2];
				Vector3d displacement = -linearDelta * TriangleNorm(p0->pos, p1->pos, p2->pos) * (TriangleSquare(p0->pos, p1->pos, p2->pos) / surface);
				/*p0->Push(p0->pos + displacement);
				p1->Push(p1->pos + displacement);
				p2->Push(p2->pos + displacement);*/
				p0->AddAcceleration(displacement);
				p1->AddAcceleration(displacement);
				p2->AddAcceleration(displacement);
			}
		}


/*		float scaleFactor = exp(1.0 / 3.0 * log(1e-3 / fabs(volume)));
		for(int i = 0; i < particlesCount; i++)
		{
			Vector3 goalPos = realMassCenter + (particles[i]->pos - realMassCenter) * scaleFactor;
			particles[i]->Push(particles[i]->pos + (goalPos - particles[i]->pos) * 1.0);
		}*/

		/*double maxAllowedDist = 0.3;
		double maxDist = 0.0;
		for(int i = 0; i < particlesCount; i++)
		{
			if((particles[i]->pos - realMassCenter).QuadLen() > sqr(maxDist))
			{
				maxDist = (particles[i]->pos - realMassCenter).Len();
			}
		}
		scaleFactor = maxAllowedDist / maxDist;
		if(scaleFactor < 1.0)
			
		{
			for(int i = 0; i < particlesCount; i++)
			{
				Vector3 goalPos = realMassCenter + (particles[i]->pos - realMassCenter) * scaleFactor;
				particles[i]->Push(particles[i]->pos + (goalPos - particles[i]->pos));
			}
		}*/


//		float scaleFactor = exp(1.0 / 3.0 * log(1.0 / (S * deformation).GetDeterminant() ));

/*		for(int i = 0; i < particlesCount; i++)
		{
			Vector3 goalPos = realMassCenter + (particles[i]->pos - realMassCenter) * scaleFactor;
			particles[i]->pos = particles[i]->pos + (goalPos - particles[i]->pos) * 0.1;
		}*/
		//resultMatrix = resultMatrix + (resultMatrix * (exp(-1.0 / 3.0 * log(fabs(resultMatrix.GetDeterminant())))) - resultMatrix) * 0.1;


		/*double scale = 1.0 / resultMatrix.GetQuadNorm();
		if(scale < 3.0)
		{
			resultMatrix = resultMatrix * scale + R * (1.0 - scale);
		}*/
		/*double maxdeform = 1.5;
		if((resultMatrix - R).GetQuadNorm() > (maxdeform))
			resultMatrix = R + (resultMatrix - R) * (maxdeform / (((resultMatrix - R).GetQuadNorm())));*/

		/*double maxdeform = 0.5;
		if((deformation - identityMatrix).GetQuadNorm() > (maxdeform))
			deformation = identityMatrix + (deformation - identityMatrix) * (maxdeform / (((deformation - identityMatrix).GetQuadNorm())));*/

		tmpx = Vector3d(S.data[0][0], S.data[1][0], S.data[2][0]);
		tmpy = Vector3d(S.data[0][1], S.data[1][1], S.data[2][1]);
		tmpz = Vector3d(S.data[0][2], S.data[1][2], S.data[2][2]);

		double mag = (tmpx ^ tmpy) * tmpz;

		Matrix3x3d invResult = resultMatrix.GetInverted();
		for(int i = 0; i < particlesCount; i++)
		{
			Vector3d goalPos = realMassCenter + resultMatrix * (destPos[i] - virtualMassCenter);
			Vector3d resPos = particles[i]->pos + (goalPos - particles[i]->pos) * stiffness;
 
			Vector3d delta = goalPos - particles[i]->pos;

			/*if(delta.SquareLen() > sqr(maxDisplacement))
			{
				Vector3d newOrigin = goalPos - delta.GetNorm() * maxDisplacement * restitution * timeStepMultiplier; 
				if(fabs(resultMatrix.GetDeterminant()) > 0.5f)
				{
					destPos[i] = invResult * (newOrigin - realMassCenter) + virtualMassCenter;
				}
			}*/

			particles[i]->AddAcceleration((goalPos - particles[i]->pos) * stiffness);
		}

    Vector3d massCenterVelocity = zeroVector3d;
    Vector3d massCenterAngularVelocity = zeroVector3d;
    for(int i = 0; i < particlesCount; i++)
    {
      massCenterVelocity += particles[i]->velocity;
      massCenterAngularVelocity += (particles[i]->pos - realMassCenter) ^ particles[i]->velocity;
    }
    massCenterVelocity /= particlesCount;
    massCenterAngularVelocity /= particlesCount;

    for(int i = 0; i < particlesCount; i++)
    {
      Vector3d desiredVelocity = massCenterVelocity;
      desiredVelocity += massCenterAngularVelocity ^ (particles[i]->pos - realMassCenter);
      particles[i]->AddAcceleration((desiredVelocity - particles[i]->velocity) * viscosity);
    }

		/*Aqq = Matrix3x3(zeroVector3, zeroVector3, zeroVector3);
		for(int i = 0; i < particlesCount; i++)
		{
			q = deformation * (destPos[i] - virtualMassCenter);
			tmp.data[0][0] = q.x * q.x; tmp.data[0][1] = q.x * q.y; tmp.data[0][2] = q.x * q.z;
			tmp.data[1][0] = q.y * q.x; tmp.data[1][1] = q.y * q.y; tmp.data[1][2] = q.y * q.z;
			tmp.data[2][0] = q.z * q.x; tmp.data[2][1] = q.z * q.y; tmp.data[2][2] = q.z * q.z;
			Aqq += tmp;
		}

		Aqq = Aqq.GetInverted();*/
	}

	static const int particlesCount = _particlesCount;
	static const int trianglesCount = _trianglesCount;

	Particle *particles[particlesCount];
	Triangle triangles[trianglesCount];

	Vector3d destPos[particlesCount];
	Matrix3x3d Aqq/*, deformation, warmStart*/;
	Matrix3x3d rotation, scale;
	double initialVolume;
	bool active;
	double stiffness;
	double compressibility;
	double restitution;
	double maxDisplacement;
  double viscosity;
};

struct Tetrahedron : public MeshlessConstraint<4, 4>
{
	Tetrahedron(){}
	Tetrahedron(Particle *_p0, Particle *_p1, Particle *_p2, Particle *_p3)
	{
		Particle *vBuf[4] = {_p0, _p1, _p2, _p3};
		Triangle tBuf[4] = {Triangle(0, 1, 2), Triangle(0, 1, 3), Triangle(1, 2, 3), Triangle(2, 0, 3)};
		Create(vBuf, tBuf);
	}
};

struct Cube : public MeshlessConstraint<8, 12>
{
	Cube(){}
	Cube(Particle *_p0, Particle *_p1, Particle *_p2, Particle *_p3,
		 Particle *_p4, Particle *_p5, Particle *_p6, Particle *_p7)
	{
		Particle *vBuf[8] = {_p0, _p1, _p2, _p3, _p4, _p5, _p6, _p7};
		Triangle tBuf[12] = {
								Triangle(0, 2, 1),
								Triangle(0, 3, 2),

								Triangle(4, 5, 6),
								Triangle(4, 6, 7),

								Triangle(1, 6, 5),
								Triangle(1, 2, 6),

								Triangle(3, 0, 4),
								Triangle(3, 4, 7),

								Triangle(2, 3, 6),
								Triangle(3, 7, 6),

								Triangle(0, 1, 5),
								Triangle(0, 5, 4)
							};
  		Create(vBuf, tBuf);
	}
};

template<int dim>
struct CubicDomain : public MeshlessConstraint<dim * dim * dim, (dim - 1) * (dim - 1) * 2 * 6>
{
	CubicDomain(){}
	CubicDomain(Particle **p)
	{
		Particle **vBuf = p;//{_p0, _p1, _p2, _p3, _p4, _p5, _p6, _p7};
		Triangle tBuf[(dim - 1) * (dim - 1) * 2 * 6];
		int trianglesCount = 0;
		for(int x = 0; x < dim - 1; x++)
			for(int z = 0; z < dim - 1; z++)
			{
				//bottom
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (z + 0) * dim + 0 * dim * dim,
													(x + 0) + (z + 1) * dim + 0 * dim * dim,
													(x + 1) + (z + 1) * dim + 0 * dim * dim);
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (z + 0) * dim + 0 * dim * dim,
													(x + 1) + (z + 1) * dim + 0 * dim * dim,
													(x + 1) + (z + 0) * dim + 0 * dim * dim);
				//top
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (z + 0) * dim + (dim - 1) * dim * dim,
													(x + 1) + (z + 1) * dim + (dim - 1) * dim * dim,
													(x + 0) + (z + 1) * dim + (dim - 1) * dim * dim);
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (z + 0) * dim + (dim - 1)* dim * dim,
													(x + 1) + (z + 0) * dim + (dim - 1) * dim * dim,
													(x + 1) + (z + 1) * dim + (dim - 1) * dim * dim);
			}

		for(int x = 0; x < dim - 1; x++)
			for(int y = 0; y < dim - 1; y++)
			{
				//front
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (0 + 0) * dim + (y + 0) * dim * dim,
													(x + 1) + (0 + 0) * dim + (y + 0) * dim * dim,
													(x + 1) + (0 + 0) * dim + (y + 1) * dim * dim);
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (0 + 0) * dim + (y + 0) * dim * dim,
													(x + 1) + (0 + 0) * dim + (y + 1) * dim * dim,
													(x + 0) + (0 + 0) * dim + (y + 1) * dim * dim);
				//back
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (dim - 1) * dim + (y + 0) * dim * dim,
													(x + 1) + (dim - 1) * dim + (y + 1) * dim * dim,
													(x + 1) + (dim - 1) * dim + (y + 0) * dim * dim);
				tBuf[trianglesCount++] = Triangle(	(x + 0) + (dim - 1) * dim + (y + 0) * dim * dim,
													(x + 0) + (dim - 1) * dim + (y + 1) * dim * dim,
													(x + 1) + (dim - 1) * dim + (y + 1) * dim * dim);
			}
		for(int z = 0; z < dim - 1; z++)
			for(int y = 0; y < dim - 1; y++)
			{
				//left
				tBuf[trianglesCount++] = Triangle(	(0 + 0) + (z + 1) * dim + (y + 0) * dim * dim,
													(0 + 0) + (z + 0) * dim + (y + 0) * dim * dim,
													(0 + 0) + (z + 0) * dim + (y + 1) * dim * dim);
				tBuf[trianglesCount++] = Triangle(	(0 + 0) + (z + 1) * dim + (y + 1) * dim * dim,
													(0 + 0) + (z + 1) * dim + (y + 0) * dim * dim,
													(0 + 0) + (z + 0) * dim + (y + 1) * dim * dim);
				//right
				tBuf[trianglesCount++] = Triangle(	(dim - 1) + (z + 1) * dim + (y + 0) * dim * dim,
													(dim - 1) + (z + 0) * dim + (y + 1) * dim * dim,
													(dim - 1) + (z + 0) * dim + (y + 0) * dim * dim);
				tBuf[trianglesCount++] = Triangle(	(dim - 1) + (z + 1) * dim + (y + 1) * dim * dim,
													(dim - 1) + (z + 0) * dim + (y + 1) * dim * dim,
													(dim - 1) + (z + 1) * dim + (y + 0) * dim * dim);
			}
		Create(vBuf, tBuf);
	}
};

typedef CubicDomain<3> Domain3x3;
typedef CubicDomain<2> Domain2x2;

class ParticleDifferentialSystem : public DifferentialSystem
{
public:
  ParticleDifferentialSystem(void (*UpdateCallback)(double time, ParticleDifferentialSystem *sys))
  {
    this->UpdateCallback = UpdateCallback;
  }

  virtual int   GetDimentionsCount()
  {
    return int(particles.size() * 6);
  }

  virtual void  GetDerivatives(double time, const double *coordinates, double *derivatives)
  {
    SetCurrCoords(coordinates);
    for(size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      particles[particleIndex].acceleration = zeroVector3d;
    }

    for(size_t cubeIndex = 0; cubeIndex < cubes.size(); cubeIndex++)
    {
      cubes[cubeIndex].Satisfy();
    }
    UpdateCallback(time, this);
    for(size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      derivatives[particleIndex * 6 + 0] = particles[particleIndex].velocity.x;
      derivatives[particleIndex * 6 + 1] = particles[particleIndex].velocity.y;
      derivatives[particleIndex * 6 + 2] = particles[particleIndex].velocity.z;

      derivatives[particleIndex * 6 + 3 + 0] = particles[particleIndex].acceleration.x;
      derivatives[particleIndex * 6 + 3 + 1] = particles[particleIndex].acceleration.y;
      derivatives[particleIndex * 6 + 3 + 2] = particles[particleIndex].acceleration.z;
    }
  }

  virtual void  GetCurrCoords(double *currCoords)
  {
    for(size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      currCoords[particleIndex * 6 + 0] = particles[particleIndex].pos.x;
      currCoords[particleIndex * 6 + 1] = particles[particleIndex].pos.y;
      currCoords[particleIndex * 6 + 2] = particles[particleIndex].pos.z;

      currCoords[particleIndex * 6 + 3 + 0] = particles[particleIndex].velocity.x;
      currCoords[particleIndex * 6 + 3 + 1] = particles[particleIndex].velocity.y;
      currCoords[particleIndex * 6 + 3 + 2] = particles[particleIndex].velocity.z;
    }
  }

  virtual void  SetCurrCoords(const double *newCoords)
  {
    for(size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      Vector3d pos, velocity;
      pos.x = newCoords[particleIndex * 6 + 0];
      pos.y = newCoords[particleIndex * 6 + 1];
      pos.z = newCoords[particleIndex * 6 + 2];

      velocity.x = newCoords[particleIndex * 6 + 3 + 0];
      velocity.y = newCoords[particleIndex * 6 + 3 + 1];
      velocity.z = newCoords[particleIndex * 6 + 3 + 2];

      particles[particleIndex].pos = pos;
      particles[particleIndex].velocity = velocity;
    }
  }

	Particle *AddParticle(Vector3d pos, Vector3d velocity = zeroVector3d)
	{
    particles.push_back(Particle(pos, velocity));
	  return &particles.back();
	}

	Cube *AddCube(Particle *p1, Particle *p2, Particle *p3, Particle *p4, Particle *p5, Particle *p6, Particle *p7, Particle *p8)
	{
    cubes.push_back(Cube(p1, p2, p3, p4, p5, p6, p7, p8));
    return &cubes.back();
	}

  std::vector<Particle> particles;
  std::vector<Cube>     cubes;
  void (*UpdateCallback)(double time, ParticleDifferentialSystem *sys);
};





ElasticSystem<Space2>::ValueType ElasticSystem<Space2>::GetGlueRiemannSolution(
  const ValueType& interiorSolution, ValueType& exteriorSolution,
  const MediumParameters& interiorParams, const MediumParameters& exteriorParams)
{
  ValueType riemannSolution;

  Scalar k0 = ((interiorSolution.GetVelocity().x - exteriorSolution.GetVelocity().x) * exteriorParams.GetZp() + interiorSolution.GetXX() - exteriorSolution.GetXX()) /
    ((interiorParams.GetZp() + exteriorParams.GetZp()) * interiorParams.GetPSpeed());

  Scalar k1 = 0;
  if (interiorParams.mju > std::numeric_limits<Scalar>::epsilon())
    k1 = ((interiorSolution.GetVelocity().y - exteriorSolution.GetVelocity().y) * exteriorParams.GetZs() + interiorSolution.GetXY() - exteriorSolution.GetXY()) /
    ((interiorParams.GetZs() + exteriorParams.GetZs()) * interiorParams.GetSSpeed());

  riemannSolution.SetTension(Tensor(
    interiorSolution.GetXX() - k0 * (interiorParams.lambda + 2 * interiorParams.mju),
    interiorSolution.GetXY() - k1 * interiorParams.mju,
    0));

  Scalar k2 = 0;
  if (interiorParams.mju > std::numeric_limits<Scalar>::epsilon())
    k2 = 1 / interiorParams.GetZs();

  riemannSolution.SetVelocity(interiorSolution.GetVelocity() -
    Vector(k0 * interiorParams.GetPSpeed(), k2 * (interiorSolution.GetXY() - riemannSolution.GetXY()))
    );
  return riemannSolution;
}

ElasticSystem<Space3>::ValueType ElasticSystem<Space3>::GetGlueRiemannSolution(
  const ValueType& interiorSolution, ValueType& exteriorSolution,
  const MediumParameters& interiorParams, const MediumParameters& exteriorParams)
{
  ValueType riemannSolution;

  Scalar k0 = ((interiorSolution.GetVelocity().x - exteriorSolution.GetVelocity().x) * exteriorParams.GetZp() + interiorSolution.GetXX() - exteriorSolution.GetXX()) /
    ((interiorParams.GetZp() + exteriorParams.GetZp()) * interiorParams.GetPSpeed());

  Scalar k1 = 0;
  if (interiorParams.mju > std::numeric_limits<Scalar>::epsilon())
    k1 = ((interiorSolution.GetVelocity().y - exteriorSolution.GetVelocity().y) * exteriorParams.GetZs() + interiorSolution.GetXY() - exteriorSolution.GetXY()) /
    ((interiorParams.GetZs() + exteriorParams.GetZs()) * interiorParams.GetSSpeed());

  Scalar k2 = 0;
  if (interiorParams.mju > std::numeric_limits<Scalar>::epsilon())
    k2 = ((interiorSolution.GetVelocity().z - exteriorSolution.GetVelocity().z) * exteriorParams.GetZs() + interiorSolution.GetXZ() - exteriorSolution.GetXZ()) /
    ((interiorParams.GetZs() + exteriorParams.GetZs()) * interiorParams.GetSSpeed());

  riemannSolution.SetTension(Tensor(
    interiorSolution.GetXX() - k0 * (interiorParams.lambda + 2 * interiorParams.mju),
    interiorSolution.GetXY() - k1 * interiorParams.mju,
    interiorSolution.GetXZ() - k2 * interiorParams.mju,
    0, 0, 0));

  Scalar k3 = 0;
  if (interiorParams.mju > std::numeric_limits<Scalar>::epsilon())
    k3 = 1 / interiorParams.GetZs();

  riemannSolution.SetVelocity(interiorSolution.GetVelocity() -
    Vector(k0 * interiorParams.GetPSpeed(), 
           k3 * (interiorSolution.GetXY() - riemannSolution.GetXY()), 
           k3 * (interiorSolution.GetXZ() - riemannSolution.GetXZ()))
    );

  return riemannSolution;
}

ElasticSystem<Space2>::ValueType ElasticSystem<Space2>::
  GetRiemannSolution(const ValueType& interiorSolution, ValueType& exteriorSolution,
    const MediumParameters& interiorParams, MediumParameters& exteriorParams,
    IndexType boundaryType, IndexType dynamicContactType)
{
  ValueType riemannSolution;

  enum {
    Contact, Boundary
  } interactionType;

  if (exteriorParams.IsZero())
  {
    interactionType = Boundary;
  } else {
    // probably contact
    interactionType = Contact;
    riemannSolution = GetGlueRiemannSolution(interiorSolution, exteriorSolution, interiorParams, exteriorParams);
  }

  if (interactionType == Boundary || riemannSolution.GetXX() > 0)
  {
    interactionType = Boundary;
      
    Eigen::Matrix<Scalar, 1, dimsCount> boundaryMatrix;
    BuildBoundaryMatrix(boundaryType, boundaryMatrix);
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      exteriorSolution.values[valueIndex] = interiorSolution.values[valueIndex] * boundaryMatrix(valueIndex);
    }
    riemannSolution = GetGlueRiemannSolution(interiorSolution, exteriorSolution, interiorParams, interiorParams);
  } else
  {
    IndexType contactType = contactDescriptions[dynamicContactType].type;
    switch (contactType)
    {
      case ContactConditions::Glue:
      {
        // do nothing
      } break;
      case ContactConditions::Glide:
      {
        riemannSolution.values[2] = 0;
      } break;
      case ContactConditions::Friction:
      {
        Scalar frictionCoeff;
        GetFrictionContactInfo(dynamicContactType, frictionCoeff);
        Scalar sign = Sgn(riemannSolution.GetXY());
        riemannSolution.values[2] = sign * std::min(fabs(riemannSolution.GetXY()), -frictionCoeff * riemannSolution.GetXX());
      } break;
      default:
        assert(0);
      break;
    }

    Scalar k2 = 0;
    if (interiorParams.mju > std::numeric_limits<Scalar>::epsilon())
      k2 = 1 / interiorParams.GetZs();

    // y-velocity correction due to sigma-xy has been changed
    riemannSolution.values[4] = interiorSolution.GetVelocity().y - k2 * (interiorSolution.GetXY() - riemannSolution.GetXY());
  }

  return riemannSolution;
}

ElasticSystem<Space3>::ValueType ElasticSystem<Space3>::
  GetRiemannSolution(const ValueType& interiorSolution, ValueType& exteriorSolution,
                     const MediumParameters& interiorParams, MediumParameters& exteriorParams,
                     IndexType boundaryType, IndexType dynamicContactType)
{
  ValueType riemannSolution;

  enum
  {
    Contact, Boundary
  } interactionType;

  if (exteriorParams.IsZero())
  {
    interactionType = Boundary;
  } else
  {
    // probably contact
    interactionType = Contact;
    riemannSolution = GetGlueRiemannSolution(interiorSolution, exteriorSolution, interiorParams, exteriorParams);
  }

  if (interactionType == Boundary || riemannSolution.GetXX() > 0)
  {
    interactionType = Boundary;

    Eigen::Matrix<Scalar, 1, dimsCount> boundaryMatrix;
    BuildBoundaryMatrix(boundaryType, boundaryMatrix);
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      exteriorSolution.values[valueIndex] = interiorSolution.values[valueIndex] * boundaryMatrix(valueIndex);
    }

    riemannSolution = GetGlueRiemannSolution(interiorSolution, exteriorSolution, interiorParams, interiorParams);
  } else
  {
    IndexType contactType = contactDescriptions[dynamicContactType].type;
    switch (contactType)
    {
    case ContactConditions::Glue:
    {
      // do nothing
    } break;
    case ContactConditions::Glide:
    {
      riemannSolution.values[3] =    // xy
      riemannSolution.values[5] = 0; // xz
    } break;
    case ContactConditions::Friction:
    {
      Scalar frictionCoeff;
      GetFrictionContactInfo(dynamicContactType, frictionCoeff);

      Scalar shearStress = sqrt(Sqr(riemannSolution.GetXY()) + Sqr(riemannSolution.GetXZ()));
      if (shearStress > -frictionCoeff * riemannSolution.GetXX())
      {
        Scalar mult = -frictionCoeff * riemannSolution.GetXX() / shearStress;
        riemannSolution.values[3] *= mult; // xy
        riemannSolution.values[5] *= mult; // xz
      }
    } break;
    default:
      assert(0);
      break;
    }

    Scalar k3 = 0;
    if (interiorParams.mju > std::numeric_limits<Scalar>::epsilon())
      k3 = 1 / interiorParams.GetZs();

    // y,z-velocities correction due to sigma-xy,xz has been changed
    riemannSolution.values[7] = interiorSolution.GetVelocity().y - k3 * (interiorSolution.GetXY() - riemannSolution.GetXY());
    riemannSolution.values[8] = interiorSolution.GetVelocity().z - k3 * (interiorSolution.GetXZ() - riemannSolution.GetXZ());
  }

  return riemannSolution;
}

void ElasticSystem<Space2>::ValueType::SetTension(Scalar xx, Scalar yy, Scalar xy)
{
  values[0] = xx;
  values[1] = yy;
  values[2] = xy;
}

void ElasticSystem<Space3>::ValueType::SetTension(Scalar xx, Scalar yy, Scalar zz, Scalar xy, Scalar yz, Scalar xz)
{
  values[0] = xx;
  values[1] = yy;
  values[2] = zz;
  values[3] = xy;
  values[4] = yz;
  values[5] = xz;
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetPressure() const
{
  return -(GetXX() + GetYY() + GetZZ()) / Scalar(3.0);
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetDeviatorSquare() const
{
  // wiki http://en.wikipedia.org/wiki/Cauchy_stress_tensor
  return Scalar(1.0 / 3) * (Sqr(GetXX() - GetYY()) + Sqr(GetXX() - GetZZ()) + Sqr(GetYY() - GetZZ())) + 
          Scalar(2.0) * (Sqr(GetXY()) + Sqr(GetXZ()) + Sqr(GetYZ()));
}

Space2::Scalar ElasticSystem<Space2>::ValueType::GetPressure() const
{
  return -(GetXX() + GetYY()) / Scalar(2.0);
}

Space2::Scalar ElasticSystem<Space2>::ValueType::GetDeviatorSquare() const
{
  return Scalar(0.5) * (Sqr(GetXX()) + Sqr(GetYY())) - 
          GetXX() * GetYY() + 
          2.0 * Sqr(GetXY());
}

Space2::Vector ElasticSystem<Space2>::ValueType::GetForce(const Vector& normal) const
{
  return Space2::Vector(
    GetXX() * normal.x + GetXY() * normal.y, 
    GetXY() * normal.x + GetYY() * normal.y);
}

Space3::Vector ElasticSystem<Space3>::ValueType::GetForce(const Vector& normal) const
{
  return Space3::Vector(
    GetXX() * normal.x + GetXY() * normal.y + GetXZ() * normal.z, 
    GetXY() * normal.x + GetYY() * normal.y + GetYZ() * normal.z,
    GetXZ() * normal.x + GetYZ() * normal.y + GetZZ() * normal.z);
}

Space2::Tensor ElasticSystem<Space2>::ValueType::GetTension() const
{
  return Space2::Tensor(GetXX(), GetXY(), GetYY());
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetZZ() const
{
  return values[2];
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetYZ() const
{
  return values[4];
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetXZ() const
{
  return values[5];
}

Space3::Tensor ElasticSystem<Space3>::ValueType::GetTension() const
{
  return Space3::Tensor(GetXX(), GetXY(), GetXZ(), GetYY(), GetYZ(), GetZZ());
}

void ElasticSystem<Space2>::BuildEdgeTransformMatrix(Vector edgeVertices[Space::NodesPerEdge], 
  MatrixXDim& transformMatrix)
{
  Vector normal  = Vector((edgeVertices[1] - edgeVertices[0]).y, -(edgeVertices[1] - edgeVertices[0]).x).GetNorm();

  transformMatrix <<
          Sqr(normal.x),         Sqr(normal.y), -Scalar(2.0) * normal.x * normal.y,         0,         0,
          Sqr(normal.y),         Sqr(normal.x),  Scalar(2.0) * normal.x * normal.y,         0,         0,
    normal.x * normal.y, -normal.x  * normal.y,      Sqr(normal.x) - Sqr(normal.y),         0,         0,
                      0,                     0,                                   0, normal.x, -normal.y,
                      0,                     0,                                   0, normal.y,  normal.x;
}

void ElasticSystem<Space2>::BuildEdgeTransformMatrixInv(Vector edgeVertices[Space::NodesPerEdge], 
  MatrixXDim& transformMatrixInv)
{
  Vector normal  = Vector((edgeVertices[1] - edgeVertices[0]).y, -(edgeVertices[1] - edgeVertices[0]).x).GetNorm();

  transformMatrixInv <<
            Sqr(normal.x),        Sqr(normal.y),  Scalar(2.0) * normal.x * normal.y,         0,        0,
            Sqr(normal.y),        Sqr(normal.x), -Scalar(2.0) * normal.x * normal.y,         0,        0,
    -normal.x  * normal.y, normal.x  * normal.y,      Sqr(normal.x) - Sqr(normal.y),         0,        0,
                        0,                    0,                                  0,  normal.x, normal.y,
                        0,                    0,                                  0, -normal.y, normal.x;
}

void ElasticSystem<Space3>::BuildFaceTransformMatrix(Vector faceVertices[Space::NodesPerFace], 
  MatrixXDim& transformMatrix)
{
  Vector normal   = ((faceVertices[1] - faceVertices[0]) ^ (faceVertices[2] - faceVertices[0])).GetNorm();
  Vector tangent0 = (faceVertices[1] - faceVertices[0]).GetNorm();
  Vector tangent1 = (normal ^ tangent0).GetNorm();

  transformMatrix <<
          Sqr(normal.x),         Sqr(tangent0.x),         Sqr(tangent1.x),           Scalar(2.0) * normal.x * tangent0.x,             Scalar(2.0) * tangent0.x * tangent1.x,           Scalar(2.0) * normal.x * tangent1.x,        0,          0,          0,
          Sqr(normal.y),         Sqr(tangent0.y),         Sqr(tangent1.y),           Scalar(2.0) * normal.y * tangent0.y,             Scalar(2.0) * tangent0.y * tangent1.y,           Scalar(2.0) * normal.y * tangent1.y,        0,          0,          0,
          Sqr(normal.z),         Sqr(tangent0.z),         Sqr(tangent1.z),           Scalar(2.0) * normal.z * tangent0.z,             Scalar(2.0) * tangent0.z * tangent1.z,           Scalar(2.0) * normal.z * tangent1.z,        0,          0,          0,
    normal.y * normal.x, tangent0.y * tangent0.x, tangent1.y * tangent1.x, normal.y * tangent0.x + normal.x * tangent0.y, tangent0.y * tangent1.x + tangent0.x * tangent1.y, normal.y * tangent1.x + normal.x * tangent1.y,        0,          0,          0,
    normal.z * normal.y, tangent0.z * tangent0.y, tangent1.z * tangent1.y, normal.z * tangent0.y + normal.y * tangent0.z, tangent0.z * tangent1.y + tangent0.y * tangent1.z, normal.z * tangent1.y + normal.y * tangent1.z,        0,          0,          0,
    normal.z * normal.x, tangent0.z * tangent0.x, tangent1.z * tangent1.x, normal.z * tangent0.x + normal.x * tangent0.z, tangent0.z * tangent1.x + tangent0.x * tangent1.z, normal.z * tangent1.x + normal.x * tangent1.z,        0,          0,          0,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.x, tangent0.x, tangent1.x,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.y, tangent0.y, tangent1.y,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.z, tangent0.z, tangent1.z;
}

void ElasticSystem<Space3>::BuildFaceTransformMatrixInv(Vector faceVertices[Space::NodesPerFace], 
  MatrixXDim& transformMatrixInv)
{
  Vector normal = ((faceVertices[1] - faceVertices[0]) ^ (faceVertices[2] - faceVertices[0])).GetNorm();
  Vector tangent0 = (faceVertices[1] - faceVertices[0]).GetNorm();
  Vector tangent1 = (normal ^ tangent0).GetNorm();
  std::swap(tangent0.x, normal.y);
  std::swap(tangent1.x, normal.z);
  std::swap(tangent1.y, tangent0.z);

  transformMatrixInv << 
          Sqr(normal.x),         Sqr(tangent0.x),         Sqr(tangent1.x),          Scalar(2.0) * normal.x  * tangent0.x,            Scalar(2.0) * tangent0.x  * tangent1.x,          Scalar(2.0) * normal.x  * tangent1.x,        0,          0,          0,
          Sqr(normal.y),         Sqr(tangent0.y),         Sqr(tangent1.y),          Scalar(2.0) * normal.y  * tangent0.y,            Scalar(2.0) * tangent0.y  * tangent1.y,          Scalar(2.0) * normal.y  * tangent1.y,        0,          0,          0,
          Sqr(normal.z),         Sqr(tangent0.z),         Sqr(tangent1.z),          Scalar(2.0) * normal.z  * tangent0.z,            Scalar(2.0) * tangent0.z  * tangent1.z,          Scalar(2.0) * normal.z  * tangent1.z,        0,          0,          0,
    normal.y * normal.x, tangent0.y * tangent0.x, tangent1.y * tangent1.x, normal.y * tangent0.x + normal.x * tangent0.y, tangent0.y * tangent1.x + tangent0.x * tangent1.y, normal.y * tangent1.x + normal.x * tangent1.y,        0,          0,          0,
    normal.z * normal.y, tangent0.z * tangent0.y, tangent1.z * tangent1.y, normal.z * tangent0.y + normal.y * tangent0.z, tangent0.z * tangent1.y + tangent0.y * tangent1.z, normal.z * tangent1.y + normal.y * tangent1.z,        0,          0,          0,
    normal.z * normal.x, tangent0.z * tangent0.x, tangent1.z * tangent1.x, normal.z * tangent0.x + normal.x * tangent0.z, tangent0.z * tangent1.x + tangent0.x * tangent1.z, normal.z * tangent1.x + normal.x * tangent1.z,        0,          0,          0,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.x, tangent0.x, tangent1.x,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.y, tangent0.y, tangent1.y, 
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.z, tangent0.z, tangent1.z;
}

void ElasticSystem<Space3>::BuildZMatrix(const MediumParameters& mediumParameters, 
  MatrixXDim& zMatrix)
{
  zMatrix <<
    mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,
                                  0, mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,
                                  0,                               0, mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                               0, -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju),
                                  0,                               0,                               0, mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                                                               0,
                                  0,                               0,                               0,                               0, mediumParameters.flowVelocity.z,                               0,                               0,           -mediumParameters.mju,                                                               0,
                                  0,                               0,                               0,                               0,                               0, mediumParameters.flowVelocity.z,           -mediumParameters.mju,                               0,                                                               0,
                                  0,                               0,                               0,                               0,                               0,        -mediumParameters.invRho, mediumParameters.flowVelocity.z,                               0,                                                               0,
                                  0,                               0,                               0,                               0,        -mediumParameters.invRho,                               0,                               0, mediumParameters.flowVelocity.z,                                                               0,
                                  0,                               0,        -mediumParameters.invRho,                               0,                               0,                               0,                               0,                               0,                                 mediumParameters.flowVelocity.z;
}

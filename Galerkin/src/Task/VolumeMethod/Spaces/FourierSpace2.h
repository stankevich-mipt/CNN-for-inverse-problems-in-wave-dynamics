template<typename GeomSpace, int Order>
struct FourierSpace
{
  typedef typename GeomSpace::IndexType IndexType;
  typedef typename GeomSpace::Scalar    Scalar;
  typedef typename GeomSpace::Vector2   Vector2;
  typedef typename GeomSpace::Vector2i  Vector2i;

  const static IndexType functionsCount = (Order + 1) * (Order + 1);
  Scalar pi;
  const static IndexType order = Order;

  FourierSpace()
  {
    pi = Scalar(3.141592653589397);
  }

  Scalar GetBasisFunctionValue(Vector2 point, IndexType functionIndex) const
  {
    Vector2i freq = GetFrequencies(functionIndex);
    return cos(pi * point.x * Scalar(freq.x) / Scalar(Order)) * cos(pi * point.y * Scalar(freq.y) / Scalar(Order));
  }

  Scalar GetBasisFunctionMaxValue(IndexType functionIndex)
  {
    return 1.0;
  }

  Scalar GetBasisFunctionDerivative(Vector2 point, Vector2i derivatives, IndexType functionIndex)
  {
    Vector2i freq = GetFrequencies(functionIndex);
    if(derivatives.x == 1 && derivatives.y == 0)
    {
      return -(pi * Scalar(freq.x) / Scalar(Order)) * sin(pi * point.x * Scalar(freq.x) / Scalar(Order)) * cos(pi * point.y * Scalar(freq.y) / Scalar(Order));
    }
    if(derivatives.x == 0 && derivatives.y == 1)
    {
      return -(pi * Scalar(freq.y) / Scalar(Order)) * cos(pi * point.x * Scalar(freq.x) / Scalar(Order)) * sin(pi * point.y * Scalar(freq.y) / Scalar(Order));
    }
    assert(0);
    return 0;
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ô1, Ô2>
  {
    Scalar res = 0;

    Scalar step = Scalar(1e-3);
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        if(x + y < Scalar(1.0))
        {
          res += GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                 GetBasisFunctionValue(Vector2(x, y), functionIndex1) *
                 step * step;
        }
      }
    }
    return res;
  }

  Vector2 ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<dÔ1/dx, Ô2>, <dÔ1/dy, Ô2>, <dÔ1/dz, Ô2>}
  {
    Vector2 res(0, 0);

    Scalar step = Scalar(1e-3);
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        if(x + y < Scalar(1.0))
        {
          res.x +=  GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                    GetBasisFunctionDerivative(Vector2(x, y), Vector2i(1, 0), functionIndex1) *
                    step * step;
          res.y +=  GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                    GetBasisFunctionDerivative(Vector2(x, y), Vector2i(0, 1), functionIndex1) *
                    step * step;
        }
      }
    }
    return res;
  }


public:
  Vector2i GetFrequencies(IndexType functionIndex)
  {
    return Vector2i(functionIndex / (Order + 1), functionIndex % (Order + 1));
  }
};

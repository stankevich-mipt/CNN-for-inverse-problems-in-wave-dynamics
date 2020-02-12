#pragma once

#include "BaseSpace.h"
#include <math.h>

template<typename Space, int Order>
struct DubinerSpace;

template<int Order>
struct DubinerSpace<Space2, Order>: public BaseSpace<Space2, Order>
{
  SPACE2_TYPEDEFS

  using BaseSpace<Space2, Order>::order;
  using BaseSpace<Space2, Order>::functionsCount;

  static const int MaxFunctionsCount = 28;

  Scalar GetBasisFunctionValue(const Vector& point, IndexType functionIndex) const
  {
    Scalar x = point.x;
    Scalar y = point.y;

    switch (functionIndex)
    {
      case 0: return 1;
      case 1: return -1 + 2 * x + y;
      case 2: return -1 + 3 * y;
      case 3: return 6 * pow(x, 2) + 6 * x * (-1 + y) + pow(-1 + y, 2);
      case 4: return (-1 + 2 * x + y) * (-1 + 5 * y); 
      case 5: return 1 - 8 * y + 10 * pow(y, 2);
      case 6: return 20 * pow(x, 3) + 30 * pow(x, 2) * (-1 + y) + 12 * x * pow(-1 + y, 2) + pow(-1 + y, 3);
      case 7: return (6 * pow(x, 2) + 6 * x * (-1 + y) + pow(-1 + y, 2)) * (-1 + 7 * y);
      case 8: return (-1 + 2 * x + y) * (1 - 12 * y + 21 * pow(y, 2));
      case 9: return -1 + 15 * y - 45 * pow(y, 2) + 35 * pow(y, 3);
      case 10: return 70 * pow(x, 4) + 140 * pow(x, 3) * (-1 + y) + 90 * pow(x, 2) * pow(-1 + y, 2) + 20 * x * pow(-1 + y, 3) + pow(-1 + y, 4);
      case 11: return (20 * pow(x, 3) + 30 * pow(x, 2) * (-1 + y) + 12 * x * pow(-1 + y, 2) + pow(-1 + y, 3)) * (-1 + 9 * y);
      case 12: return (6 * pow(x, 2) + 6 * x * (-1 + y) + pow(-1 + y, 2)) * (1 - 16 * y + 36 * pow(y, 2));
      case 13: return (-1 + 2 * x + y) * (-1 + 21 * y - 84 * pow(y, 2) + 84 * pow(y, 3));
      case 14: return 1 - 24 * y + 126 * pow(y, 2) - 224 * pow(y, 3) + 126 * pow(y, 4);
      case 15: return 252 * pow(x, 5) + 630 * pow(x, 4) * (-1 + y) + 560 * pow(x, 3) * pow(-1 + y, 2) + 210 * pow(x, 2) * pow(-1 + y, 3) + 30 * x * pow(-1 + y, 4) + pow(-1 + y, 5);
      case 16: return (70 * pow(x, 4) + 140 * pow(x, 3) * (-1 + y) + 90 * pow(x, 2) * pow(-1 + y, 2) + 20 * x * pow(-1 + y, 3) + pow(-1 + y, 4)) * (-1 + 11 * y);
      case 17: return (20 * pow(x, 3) + 30 * pow(x, 2) * (-1 + y) + 12 * x * pow(-1 + y, 2) + pow(-1 + y, 3)) * (1 - 20 * y + 55 * pow(y, 2));
      case 18: return (6 * pow(x, 2) + 6 * x * (-1 + y) + pow(-1 + y, 2)) * (-1 + 27 * y - 135 * pow(y, 2) + 165 * pow(y, 3));
      case 19: return (-1 + 2 * x + y) * (1 - 32 * y + 216 * pow(y, 2) - 480 * pow(y, 3) + 330 * pow(y, 4));
      case 20: return -1 + 35 * y - 280 * pow(y, 2) + 840 * pow(y, 3) - 1050 * pow(y, 4) + 462 * pow(y, 5);
      case 21: return 924 * pow(x, 6) + 2772 * pow(x, 5) * (-1 + y) + 3150 * pow(x, 4) * pow(-1 + y, 2) + 1680 * pow(x, 3) * pow(-1 + y, 3) + 420 * pow(x, 2) * pow(-1 + y, 4) + 42 * x * pow(-1 + y, 5) + pow(-1 + y, 6);
      case 22: return (252 * pow(x, 5) + 630 * pow(x, 4) * (-1 + y) + 560 * pow(x, 3) * pow(-1 + y, 2) + 210 * pow(x, 2) * pow(-1 + y, 3) + 30 * x * pow(-1 + y, 4) + pow(-1 + y, 5)) * (-1 + 13 * y);
      case 23: return (70 * pow(x, 4) + 140 * pow(x, 3) * (-1 + y) + 90 * pow(x, 2) * pow(-1 + y, 2) + 20 * x * pow(-1 + y, 3) + pow(-1 + y, 4)) * (1 - 24 * y + 78 * pow(y, 2));
      case 24: return (20 * pow(x, 3) + 30 * pow(x, 2) * (-1 + y) + 12 * x * pow(-1 + y, 2) + pow(-1 + y, 3)) * (-1 + 33 * y - 198 * pow(y, 2) + 286 * pow(y, 3));
      case 25: return (6 * pow(x, 2) + 6 * x * (-1 + y) + pow(-1 + y, 2)) * (1 - 40 * y + 330 * pow(y, 2) - 880 * pow(y, 3) + 715 * pow(y, 4));
      case 26: return (-1 + 2 * x + y) * (-1 + 45 * y - 450 * pow(y, 2) + 1650 * pow(y, 3) - 2475 * pow(y, 4) + 1287 * pow(y, 5));
      case 27: return 1 - 48 * y + 540 * pow(y, 2) - 2400 * pow(y, 3) + 4950 * pow(y, 4) - 4752 * pow(y, 5) + 1716 * pow(y, 6);
      default:
        assert(0);
      break;
    }
    return 0;
  }
};

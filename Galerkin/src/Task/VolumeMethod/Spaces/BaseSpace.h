#pragma once

#include "../../../Maths/Spaces.h"

template <typename Space, int Order>
struct BaseSpace;

template <int Order>
struct BaseSpace<Space2, Order>
{
  SPACE2_TYPEDEFS
  const static IndexType functionsCount = (Order + 1) * (Order + 2) / 2;
  const static IndexType order = Order;
};

template <int Order>
struct BaseSpace<Space3, Order>
{
  SPACE3_TYPEDEFS
  const static IndexType functionsCount = (Order + 1) * (Order + 2) * (Order + 3) / 6;
  const static IndexType order = Order;
};

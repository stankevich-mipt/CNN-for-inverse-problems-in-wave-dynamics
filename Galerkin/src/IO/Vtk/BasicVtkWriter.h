#pragma once

#include <string>
#include <fstream>
#include <algorithm>

#include "../../Maths/Spaces.h"

template <typename Space>
class BasicVtkWriter
{
  SPACE_TYPEDEFS

  protected:
    static std::string TypeName();
    static void WriteNumber(std::fstream&, Scalar); // ascii format
    void WriteTuple(std::fstream&,
                    Scalar, Scalar);
    void WriteTuple(std::fstream&,
                    Scalar, Scalar, Scalar);
    void WriteTuple(std::fstream&,
                    Scalar, Scalar, Scalar, Scalar);
};

template <typename Space>
void BasicVtkWriter<Space>::WriteNumber(std::fstream& file, Scalar value)
{
  // binary
  swap_endian(value);
  file.write(static_cast<const char *>(&value), sizeof(Scalar));
  
  // ascii
  // file << value << " ";
}

template <typename Space>
void BasicVtkWriter<Space>::WriteTuple(std::fstream& file,
                                       Scalar x, Scalar y)
{
  WriteNumber(file, x);
  WriteNumber(file, y);
}

template <typename Space>
void BasicVtkWriter<Space>::WriteTuple(std::fstream& file,
  Scalar x, Scalar y, Scalar z)
{
  WriteNumber(file, x);
  WriteNumber(file, y);
  WriteNumber(file, z);
}

template <typename Space>
void BasicVtkWriter<Space>::WriteTuple(std::fstream& file,
  Scalar x, Scalar y, Scalar z, Scalar w)
{
  WriteNumber(file, x);
  WriteNumber(file, y);
  WriteNumber(file, z);
  WriteNumber(file, w);
}

template <typename Space>
std::string BasicVtkWriter<Space>::TypeName()
{
  if (sizeof(Scalar) == sizeof(float)) 
  {
    return "float";
  } else 
  if (sizeof(Scalar) == sizeof(double)) 
  {
    return "double"; 
  } else 
  if (sizeof(Scalar) == sizeof(long double)) 
  {
    return "long double";
  } else 
  {
    return "unknown type";
  }
}

template<typename Space>
void ComputeSnapshotAABB(
  const typename Space::Vector& origin, const typename Space::Vector& spacing, 
  typename Space::AABB aabb, typename Space::IntAABB* iaabb)
{
  typedef typename Space::Scalar Scalar;

  Scalar eps = std::numeric_limits<Scalar>::epsilon();

  for (typename Space::IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
  {
    assert(spacing.Get(dimIndex) > Scalar(0.0));
    Scalar min = (aabb.boxPoint1[dimIndex] - origin.Get(dimIndex)) / spacing.Get(dimIndex);
    iaabb->boxPoint1[dimIndex] = static_cast<int>(min - fabs(min * eps));

    // no more than 2 iterations
    while (origin.Get(dimIndex) + spacing.Get(dimIndex) * Scalar(iaabb->boxPoint1.Get(dimIndex)) + eps < aabb.boxPoint1.Get(dimIndex))
    {
      ++iaabb->boxPoint1[dimIndex];
    }

    typename Space::Scalar max = (aabb.boxPoint2[dimIndex] - origin.Get(dimIndex)) / spacing.Get(dimIndex);
    iaabb->boxPoint2[dimIndex] = static_cast<int>(max + fabs(max * eps));
    while (origin.Get(dimIndex) + spacing.Get(dimIndex) * Scalar(iaabb->boxPoint2.Get(dimIndex)) > aabb.boxPoint2.Get(dimIndex) + eps)
    {
      --iaabb->boxPoint2[dimIndex];
    }
  }
}

template<typename Space>
typename Space::IndexType ComputeSnapshotSize(
  const typename Space::Vector& origin, const typename Space::Vector& spacing, 
  const typename Space::AABB aabb)
{
  typename Space::IntAABB iaabb;
  ComputeSnapshotAABB<Space>(origin, spacing, aabb, &iaabb);

  return ((iaabb.boxPoint2 - iaabb.boxPoint1).ComponentAbs() + Space::IndexVector::one()).GetVolume();
}

template<typename Space>
typename Space::IndexType ComputeSnapshotSize(
  const typename Space::Vector& origin, const typename Space::Vector& spacing, 
  const typename Space::AABB aabb, typename Space::IntAABB* iaabb)
{
  ComputeSnapshotAABB<Space>(origin, spacing, aabb, iaabb);
  return ((iaabb->boxPoint2 - iaabb->boxPoint1).ComponentAbs() + Space::IndexVector::one()).GetVolume();
}

template<typename Space>
typename Space::IndexVector ComputeSnapshotAABBSize(
  const typename Space::Vector& origin, const typename Space::Vector& spacing, 
  const typename Space::AABB aabb)
{
  typename Space::IntAABB iaabb;
  ComputeSnapshotAABB<Space>(origin, spacing, aabb, &iaabb);
  return (iaabb.boxPoint2 - iaabb.boxPoint1).ComponentAbs() + Space::IndexVector::one();
}

template<typename Space>
typename Space::IndexVector ComputeSnapshotAABBSize(
  const typename Space::Vector& origin, const typename Space::Vector& spacing, 
  const typename Space::AABB aabb, typename Space::IntAABB* iaabb)
{
  ComputeSnapshotAABB<Space>(origin, spacing, aabb, iaabb);
  return (iaabb->boxPoint2 - iaabb->boxPoint1).ComponentAbs() + Space::IndexVector::one();
}

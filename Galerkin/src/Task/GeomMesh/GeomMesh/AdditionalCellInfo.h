#pragma once

#include "../../../Maths/Spaces.h"

template<typename Space>
struct AdditionalCellInfoBase;

template<>
struct AdditionalCellInfoBase<Space2>
{
  SPACE2_TYPEDEFS
  struct NeighbourEdge
  {
    NeighbourEdge(): correspondingCellIndex(IndexType(-1)),
                     correspondingEdgeNumber(IndexType(-1)),
                     // -1 - indefinite value
                     interactionType(IndexType(-1))
    {}
    IndexType correspondingCellIndex;
    IndexType correspondingEdgeNumber;
    IndexType interactionType;
  };
  NeighbourEdge neighbouringEdges[Space2::EdgesPerCell];
};

template<>
struct AdditionalCellInfoBase<Space3>
{
  SPACE3_TYPEDEFS
  struct NeighbourFace
  {
    NeighbourFace(): correspondingCellIndex(IndexType(-1)),
                     correspondingFaceNumber(IndexType(-1)),
                     // -1 - indefinite value
                     interactionType(IndexType(-1)),
                     orientation(IndexType(-1))
    {}
    IndexType correspondingCellIndex;
    IndexType correspondingFaceNumber;
    IndexType interactionType;
    IndexType orientation;
  };
  NeighbourFace neighbouringFaces[Space3::FacesPerCell];
};

template <typename Space>
struct AdditionalCellInfo: public AdditionalCellInfoBase<Space>
{
  SPACE_TYPEDEFS

  template <typename T>
  struct AuxInfo
  {
    SPACE_TYPEDEFS
    typedef T DataType;
    AuxInfo(): data(nullptr), count(0)
    {}
    DataType* data;
    IndexType count;
  };
};

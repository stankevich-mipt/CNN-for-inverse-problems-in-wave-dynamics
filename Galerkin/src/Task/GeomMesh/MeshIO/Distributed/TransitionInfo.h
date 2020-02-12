#pragma once

#include "../../GeomMesh/GeomMesh.h"
#include <vector>

template <typename Space>
struct TransitionInfo
{
  SPACE3_TYPEDEFS
  typedef typename GeomMesh<Space>::Cell Cell;
  struct TransitionNode
  {
    IndexType nativeIndex;
    IndexType targetIndex;
    friend std::ostream& operator<<(std::ostream& os, const TransitionNode& transitionNode)
    {
      os << transitionNode.nativeIndex << " " << transitionNode.targetIndex;
      return os;
    }
    friend std::istream& operator>>(std::istream& os, TransitionNode& transitionNode)
    {
      os >> transitionNode.nativeIndex >> transitionNode.targetIndex;
      return os;
    }
  };

  TransitionInfo& operator=(const TransitionInfo& other) 
  {
    nodesIndices = other.nodesIndices;
    cells = other.cells;
    transitionNodes = other.transitionNodes;
    return *this;
  }

  std::vector<IndexType>      nodesIndices;
  std::vector<Cell>           cells;
  std::vector<TransitionNode> transitionNodes;
};

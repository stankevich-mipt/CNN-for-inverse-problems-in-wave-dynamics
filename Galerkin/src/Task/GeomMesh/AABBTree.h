#pragma once
#include <vector>
#include <algorithm>
#include "../../Maths/Spaces.h"

template <typename Space, typename UserData>
class AABBTree
{
public:
  SPACE_TYPEDEFS

  AABBTree()
  {
    firstFreeNodeIndex = -1;
    rootNodeIndex = -1;
  }

  int InsertNode(AABB aabb)
  {
    if(rootNodeIndex == -1) 
    {
      rootNodeIndex = AllocateNode();
      nodes[rootNodeIndex].aabb = aabb;
      nodes[rootNodeIndex].isLeaf = 1;
      nodes[rootNodeIndex].parentIndex = -1;
      nodes[rootNodeIndex].height = 0;
      return rootNodeIndex;
    }

    int currNodeIndex = FindOptimalLocation(aabb);

    int newLeafIndex   = AllocateNode();
    int newBranchIndex = AllocateNode();

    if(currNodeIndex == rootNodeIndex)
      rootNodeIndex = newBranchIndex;

    nodes[newBranchIndex].childrenIndices[0] = currNodeIndex;
    nodes[newBranchIndex].childrenIndices[1] = newLeafIndex;
    nodes[newBranchIndex].parentIndex = nodes[currNodeIndex].parentIndex;
    nodes[newBranchIndex].isLeaf = 0;

    if(nodes[newBranchIndex].parentIndex != -1)
    {
      int parentIndex = nodes[newBranchIndex].parentIndex;
      if(nodes[parentIndex].childrenIndices[0] == currNodeIndex) 
        nodes[parentIndex].childrenIndices[0] = newBranchIndex;
      if(nodes[parentIndex].childrenIndices[1] == currNodeIndex) 
        nodes[parentIndex].childrenIndices[1] = newBranchIndex;
    }

    nodes[currNodeIndex].parentIndex = newBranchIndex;

    nodes[newLeafIndex].aabb = aabb;
    nodes[newLeafIndex].aabb.Expand(1.1f);
    nodes[newLeafIndex].parentIndex = newBranchIndex;
    nodes[newLeafIndex].isLeaf = 1;
    nodes[newLeafIndex].height = 0;

    UpdateNodeHierarchy(newBranchIndex);

    return newLeafIndex;
  }

  void SetUserData(int nodeIndex, UserData data)
  {
    nodes[nodeIndex].userData = data;
  }
  UserData GetUserData(int nodeIndex) const
  {
    return nodes[nodeIndex].userData;
  }

  void UpdateNode(int nodeIndex, AABB aabb)
  {
    assert(nodes[nodeIndex].IsLeaf());
    if(nodes[nodeIndex].aabb.Includes(aabb)) return;

    RemoveNode(nodeIndex); //nodeIndex will be placed to the top of free nodes pool
    int newIndex = InsertNode(aabb);
    assert(newIndex == nodeIndex); //ensure that nodeIndex has not changed and user data is still valid
  }

  void RemoveNode(int nodeIndex)
  {
    assert(nodes[nodeIndex].IsLeaf());
    int parentIndex = nodes[nodeIndex].parentIndex;

    int siblingIndex = -1;
    if(parentIndex != -1)
    {
      if(nodes[parentIndex].childrenIndices[0] == nodeIndex)
        siblingIndex = nodes[parentIndex].childrenIndices[1];
      else
        siblingIndex = nodes[parentIndex].childrenIndices[0];

      int grandParentIndex = nodes[parentIndex].parentIndex;
      nodes[siblingIndex].parentIndex = grandParentIndex;

      if(parentIndex == rootNodeIndex)
        rootNodeIndex = siblingIndex;

      if(grandParentIndex != -1)
      {
        if(nodes[grandParentIndex].childrenIndices[0] == parentIndex)
          nodes[grandParentIndex].childrenIndices[0] = siblingIndex;
        else
          nodes[grandParentIndex].childrenIndices[1] = siblingIndex;

        UpdateNodeHierarchy(grandParentIndex);
      }

      DeallocateNode(parentIndex);
    }
    DeallocateNode(nodeIndex);
    if(nodeIndex == rootNodeIndex)
      rootNodeIndex = -1;
  }

  template<typename Functor>
  void FindCollisions(AABB aabb, Functor &functor) const
  {
    if(rootNodeIndex != -1)
      FindCollisionsRecursive(aabb, functor, rootNodeIndex);
  }

  template<typename OtherTree, typename Functor>
  void FindTreeCollisions(const OtherTree &otherTree, Functor &functor) const
  {
    if(rootNodeIndex != -1 && otherTree.rootNodeIndex != -1)
      FindTreeCollisionsRecursive(otherTree, functor, rootNodeIndex, otherTree.rootNodeIndex);
  }

  void Erase()
  {
    rootNodeIndex = -1;
    firstFreeNodeIndex = nodes.size() > 0 ? 0 : -1;
    for(size_t nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++)
    {
      nodes[nodeIndex].nextIndex = (nodeIndex + 1) < nodes.size() ? nodeIndex + 1 : -1;
    };
  }

private:
  struct Node
  {
    AABB aabb;
    union
    {
      int parentIndex; //if in hierarchy
      int nextIndex;   //if free
    };

    union
    {
      int childrenIndices[2]; //only branches have children
      UserData userData;   //only leaves have userdata
    };
    int height;
    bool isLeaf;

    bool IsLeaf() const
    {
      return isLeaf;
    }
  };

public:
  std::vector<Node> nodes;
  int rootNodeIndex;
  int firstFreeNodeIndex;
private:

  int FindOptimalLocation(AABB aabb)
  {
    int currNodeIndex = rootNodeIndex;

    while(!nodes[currNodeIndex].IsLeaf())
    {
      /*AABB aabbs[2];
      for(int childNumber = 0; childNumber < 2; childNumber++)
      {
        int childIndex = nodes[currNodeIndex].childrenIndices[childNumber];
        aabbs[childNumber] = nodes[childIndex].aabb;
        aabbs[childNumber].Expand(aabb);
      }
      if(aabbs[0].Square() < aabbs[1].Square())
        currNodeIndex = nodes[currNodeIndex].childrenIndices[0];
      else
        currNodeIndex = nodes[currNodeIndex].childrenIndices[1];*/
      int childrenIndices[2];
      childrenIndices[0] = nodes[currNodeIndex].childrenIndices[0];
      childrenIndices[1] = nodes[currNodeIndex].childrenIndices[1];


      Scalar square = nodes[currNodeIndex].aabb.Square();

      AABB sumAABB = nodes[currNodeIndex].aabb;
      sumAABB.Expand(aabb);

      Scalar sumSquare = sumAABB.Square();


      // Cost of creating a new parent for this node and the new leaf
      Scalar baseCost = 2.0f * sumSquare;


      // Minimum cost of pushing the leaf further down the tree
      Scalar inheritanceCost = 2.0f * (sumSquare - square);


      // Cost of descending into child1
      Scalar childrenCost[2];
      for(int childNumber = 0; childNumber < 2; childNumber++)
      {
        int childIndex = childrenIndices[childNumber];
        if (nodes[childIndex].IsLeaf())
        {
          AABB sumAABB = nodes[childIndex].aabb;
          sumAABB.Expand(aabb);
          childrenCost[childNumber] = sumAABB.Square() + inheritanceCost;
        } else
        {
          AABB sumAABB = nodes[childIndex].aabb;
          sumAABB.Expand(aabb);
          Scalar oldArea = nodes[childIndex].aabb.Square();
          Scalar newArea = sumAABB.Square();
          childrenCost[childNumber] = (newArea - oldArea) + inheritanceCost;
        }
      }

      // Descend according to the minimum cost.
      if (baseCost < childrenCost[0] && baseCost < childrenCost[1])
      {
        return currNodeIndex;
      }

      // Descend
      if (childrenCost[0] < childrenCost[1])
      {
        currNodeIndex = childrenIndices[0];
      }
      else
      {
        currNodeIndex = childrenIndices[1];
      }
    }
    return currNodeIndex;
  }
  void UpdateNodeHierarchy(int nodeIndex)
  {
    int currNode = nodeIndex;
    while(currNode != -1)
    {
      nodes[currNode].height = std::max(
        nodes[nodes[currNode].childrenIndices[0]].height, 
        nodes[nodes[currNode].childrenIndices[1]].height) + 1;

      nodes[currNode].aabb =      nodes[nodes[currNode].childrenIndices[0]].aabb;
      nodes[currNode].aabb.Expand(nodes[nodes[currNode].childrenIndices[1]].aabb);

      currNode = nodes[currNode].parentIndex;
    }
  }

  template<typename Functor>
  void FindCollisionsRecursive(const AABB &aabb, Functor& functor, int currNode) const
  {
    if(!aabb.Intersects(nodes[currNode].aabb)) return;

    if(!nodes[currNode].IsLeaf())
    {
      for(int childNumber = 0; childNumber < 2; childNumber++)
        FindCollisionsRecursive(aabb, functor, nodes[currNode].childrenIndices[childNumber]);
    }else
    {
      functor(currNode);
    }
  }

  template<typename OtherTree, typename Functor>
  void FindTreeCollisionsRecursive(const OtherTree &otherTree, Functor &functor, int currNode0, int currNode1) const
  {
    if(!nodes[currNode0].aabb.Intersects(otherTree.nodes[currNode1].aabb)) return;

    if(!nodes[currNode0].IsLeaf() && !otherTree.nodes[currNode1].IsLeaf())
    {
      for(int childNumber0 = 0; childNumber0 < 2; childNumber0++)
        for(int childNumber1 = 0; childNumber1 < 2; childNumber1++)
          FindTreeCollisionsRecursive(otherTree, functor, 
            nodes[currNode0].childrenIndices[childNumber0], 
            otherTree.nodes[currNode1].childrenIndices[childNumber1]);
    }else
    if(nodes[currNode0].IsLeaf() && !otherTree.nodes[currNode1].IsLeaf())
    {
      for(int childNumber1 = 0; childNumber1 < 2; childNumber1++)
        FindTreeCollisionsRecursive(otherTree, functor, currNode0, otherTree.nodes[currNode1].childrenIndices[childNumber1]);
    }else
    if(!nodes[currNode0].IsLeaf() && otherTree.nodes[currNode1].IsLeaf())
    {
      for(int childNumber0 = 0; childNumber0 < 2; childNumber0++)
        FindTreeCollisionsRecursive(otherTree, functor, nodes[currNode0].childrenIndices[childNumber0], currNode1);
    }else
    {
      functor(currNode0, currNode1);
    }
  }

  void DeallocateNode(int nodeIndex)
  {
    nodes[nodeIndex].height = -1;
    nodes[nodeIndex].nextIndex = firstFreeNodeIndex;
    firstFreeNodeIndex = nodeIndex;
  }
  int AllocateNode()
  {
    if(firstFreeNodeIndex == -1)
    {
      firstFreeNodeIndex = nodes.size();
      Node newbie;
      newbie.height = -1;
      newbie.nextIndex = -1;
      nodes.push_back(newbie);
    }
    int poolTop = firstFreeNodeIndex;
    firstFreeNodeIndex = nodes[firstFreeNodeIndex].nextIndex;
    return poolTop;
  }
};

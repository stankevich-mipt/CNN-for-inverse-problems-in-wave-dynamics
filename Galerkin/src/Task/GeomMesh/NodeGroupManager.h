#pragma once

#include <vector>

template <typename Space, typename VolumeMeshType>
class NodeGroupManager
{
public:
  SPACE_TYPEDEFS

  NodeGroupManager(IndexType maxNodeGroupSize): 
    maxNodeGroupSize(maxNodeGroupSize), volumeMesh(nullptr), groupsCount(0)
  {
  }

  void SetVolumeMesh(VolumeMeshType* volumeMesh)
  {
    this->volumeMesh = volumeMesh;

    nodeGroupPool.resize(volumeMesh->nodes.size() * maxNodeGroupSize, 0);
    nodeGroupSizes.resize(volumeMesh->nodes.size(), 0);
    nodeGroupIndices.resize(volumeMesh->nodes.size(), IndexType(-1));
  }

  void BuildNodeGroups()
  {
    for (IndexType nodeIndex = 0; nodeIndex < volumeMesh->nodes.size(); ++nodeIndex)
    {
      if (nodeGroupIndices[nodeIndex] == IndexType(-1))
      {
        std::vector<IndexType> pool;
        pool.reserve(maxNodeGroupSize);

        nodeGroupSizes[groupsCount] = volumeMesh->GetNodeGroup(nodeIndex, pool);

        for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSizes[groupsCount]; ++nodeNumber)
        {
          nodeGroupIndices[pool[nodeNumber]] = groupsCount;
          nodeGroupPool[groupsCount * maxNodeGroupSize + nodeNumber] = pool[nodeNumber];
        }

        ++groupsCount;
      }
    }
  }

  void UpdateNode(IndexType nodeIndex)
  {
    IndexType groupIndex = nodeGroupIndices[nodeIndex];
    if (groupIndex == IndexType(-1)) return;

    std::vector<IndexType> group;
    group.reserve(nodeGroupSizes[groupIndex]);

    for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSizes[groupIndex]; ++nodeNumber)
    {
      group.push_back(nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber]);
    }

    RemoveGroup(groupIndex);

    for (IndexType nodeIndex : group)
    {
      if (nodeGroupIndices[nodeIndex] == IndexType(-1))
      {
        std::vector<IndexType> pool;
        pool.reserve(maxNodeGroupSize);

        nodeGroupSizes[groupsCount] = volumeMesh->GetNodeGroup(nodeIndex, pool);

        for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSizes[groupsCount]; ++nodeNumber)
        {
          nodeGroupIndices[pool[nodeNumber]] = groupsCount;
          nodeGroupPool[groupsCount * maxNodeGroupSize + nodeNumber] = pool[nodeNumber];
        }

        ++groupsCount;
      }
    }
  }

  void UpdateFace(IndexType cellIndex, IndexType faceNumber)
  {
    IndexType faceNodeIndices[Space::NodesPerFace];
    volumeMesh->GetCellFaceNodes(cellIndex, faceNumber, faceNodeIndices);

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
    {
      UpdateNode(faceNodeIndices[nodeNumber]);
    }
  }

  void RemoveNode(IndexType nodeIndex)
  {
    IndexType groupIndex = nodeGroupIndices[nodeIndex];
    assert(groupIndex != IndexType(-1));

    IndexType nodeNumber = 0;
    while (nodeNumber < nodeGroupSizes[groupIndex] && nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber] != nodeIndex)
      ++nodeNumber;

    assert(nodeNumber < nodeGroupSizes[groupIndex]);

    std::swap(nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber], nodeGroupPool[groupIndex * maxNodeGroupSize + nodeGroupSizes[groupIndex] - 1]);
    nodeGroupPool[groupIndex * maxNodeGroupSize + nodeGroupSizes[groupIndex] - 1] = 0;
    --nodeGroupSizes[groupIndex];
    nodeGroupIndices[nodeIndex] = IndexType(-1);
  }

  void RemoveCell(IndexType cellIndex)
  {
    IndexType nodeIndices[Space::NodesPerCell];
    volumeMesh->GetFixedCellIndices(cellIndex, nodeIndices);

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      RemoveNode(nodeIndices[nodeNumber]);
    }
  }

  IndexType GetGroupSize(IndexType nodeIndex) const
  {
    IndexType groupIndex = nodeGroupIndices[nodeIndex];
    // assert(groupIndex != IndexType(-1));
    return groupIndex != IndexType(-1) ? nodeGroupSizes[groupIndex] : 0;
  }

  const IndexType* GetGroup(IndexType nodeIndex) const
  {
    IndexType groupIndex = nodeGroupIndices[nodeIndex];
    return groupIndex != IndexType(-1) ? nodeGroupPool.data() + groupIndex * maxNodeGroupSize : 0;
  }

private:
  const IndexType maxNodeGroupSize;
  VolumeMeshType* volumeMesh;
  std::vector<IndexType> nodeGroupPool;
  std::vector<IndexType> nodeGroupIndices;
  std::vector<IndexType> nodeGroupSizes;
  IndexType groupsCount;

  void RemoveGroup(IndexType groupIndex)
  {
    for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSizes[groupIndex]; ++nodeNumber)
    {
      IndexType nodeIndex = nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber];
      nodeGroupIndices[nodeIndex] = IndexType(-1);
    }

    for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSizes[groupsCount - 1]; ++nodeNumber)
    {
      IndexType nodeIndex = nodeGroupPool[(groupsCount - 1) * maxNodeGroupSize + nodeNumber];
      nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber] = nodeIndex;
    }

    nodeGroupSizes[groupIndex] = nodeGroupSizes[groupsCount - 1];
    nodeGroupSizes[groupsCount - 1] = 0;
    --groupsCount;

    for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSizes[groupIndex]; ++nodeNumber)
    {
      IndexType nodeIndex = nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber];
      nodeGroupIndices[nodeIndex] = groupIndex;
    }
  }
};

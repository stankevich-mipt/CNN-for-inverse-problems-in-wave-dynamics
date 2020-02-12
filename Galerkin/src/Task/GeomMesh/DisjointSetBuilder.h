#include <vector>

template<typename IndexType>
struct DisjointSetUnion
{
  DisjointSetUnion(IndexType setsCount)
  {
    setIndices.resize(setsCount);
    for(IndexType setIndex = 0; setIndex < setsCount; ++setIndex)
    {
      setIndices[setIndex] = setIndex;
    }
  };

  IndexType GetSetIndex(IndexType elementIndex)
  {
    if(setIndices[elementIndex] == elementIndex) return elementIndex;
    setIndices[elementIndex] = GetSetIndex(setIndices[elementIndex]);
    return setIndices[elementIndex];
  }

  void MergeElements(IndexType elementIndex0, IndexType elementIndex1)
  {
    setIndices[GetSetIndex(elementIndex0)] = setIndices[GetSetIndex(elementIndex1)];
  }

  std::vector<IndexType> setIndices;
};

template<typename T>
struct Pair
{
  T values[2];
};

template<typename T, typename IndexType>
struct SetElement
{
  T value;
  IndexType pairIndex;
  static bool CompareValue(SetElement elem0, SetElement elem1)
  {
    return elem0.value < elem1.value;
  }

  static bool CompareIndex(SetElement elem0, SetElement elem1)
  {
    return elem0.pairIndex < elem1.pairIndex;
  }
};

template<typename T, typename IndexType>
struct GroupBuilder
{
  GroupBuilder(){}

  void LoadPairs(Pair<T> *pairs, IndexType pairsCount)
  {
    if(pairsCount == 0) return;

    std::vector<SetElement<T, IndexType> > elements;

    elements.resize(pairsCount * 2);
    for(IndexType pairIndex = 0; pairIndex < pairsCount; ++pairIndex)
    {
      for(IndexType valueIndex = 0; valueIndex < 2; ++valueIndex)
      {
        SetElement<T, IndexType> newbie;
        newbie.pairIndex = pairIndex;
        newbie.value = pairs[pairIndex].values[valueIndex];
        elements[pairIndex * 2 + valueIndex] = newbie;
      }
    }

    /*for(IndexType i = 0; i < elements.size(); i++)
    {
      std::cout << "(" << elements[i].value.edgeIndex << ", " << elements[i].pairIndex << ") ";
    }*/

    DisjointSetUnion<IndexType> dsu(pairsCount);
    std::sort(elements.begin(), elements.end(), SetElement<T, IndexType>::CompareValue);
    IndexType size = IndexType(elements.size());
    for(IndexType i = 1; i < size; ++i)
    {
      if(elements[i - 1].value == elements[i].value)
      {
        dsu.MergeElements(elements[i - 1].pairIndex, elements[i].pairIndex);
      }
    }


    std::vector<IndexType> frequencies;
    frequencies.resize(pairsCount);

    std::vector<IndexType> pairToGroup;
    pairToGroup.resize(pairsCount);

    for(IndexType pairIndex = 0; pairIndex < pairsCount; ++pairIndex)
    {
      ++frequencies[dsu.GetSetIndex(pairIndex)];
    }

    IndexType uniqueGroupsCount = 0;
    for(IndexType pairIndex = 0; pairIndex < pairsCount; ++pairIndex)
    {
      if(frequencies[pairIndex] != 0) 
      {
        pairToGroup[pairIndex] = uniqueGroupsCount++;
      }
    }

    pairGroupIndices.resize(pairsCount);
    groupInfos.resize(uniqueGroupsCount);
    for(IndexType pairIndex = 0; pairIndex < pairsCount; ++pairIndex)
    {
      ++groupInfos[pairToGroup[dsu.GetSetIndex(pairIndex)]].size;
    }

    IndexType test0 = groupInfos[0].size;

    IndexType currGroupOffset = 0;
    for(IndexType groupIndex = 0; groupIndex < uniqueGroupsCount; ++groupIndex)
    {
      groupInfos[groupIndex].offset = currGroupOffset;
      currGroupOffset += groupInfos[groupIndex].size;
      groupInfos[groupIndex].size = 0; //will be recalculated
    }

    for(IndexType pairIndex = 0; pairIndex < pairsCount; ++pairIndex)
    {
      IndexType groupIndex = pairToGroup[dsu.GetSetIndex(pairIndex)];
      pairGroupIndices[groupInfos[groupIndex].offset + groupInfos[groupIndex].size++] = pairIndex;
    }

  }
  IndexType GetGroupsCount()
  {
    return groupInfos.size();
  }
  IndexType GetGroupSize(IndexType groupIndex)
  {
    return groupInfos[groupIndex].size;
  }
  IndexType GetGroupElement(IndexType groupIndex, IndexType elementIndex)
  {
    return pairGroupIndices[groupInfos[groupIndex].offset + elementIndex];
  }
private:

  struct GroupInfo
  {
    IndexType size;
    IndexType offset;
  };

  std::vector<GroupInfo> groupInfos;
  std::vector<IndexType> pairGroupIndices;
};

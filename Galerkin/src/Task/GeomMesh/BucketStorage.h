#pragma once
#include <assert.h>

template <typename IndexType>
class BucketStorage
{
public:
  typedef const IndexType* ElementIterator;

  BucketStorage():
    bucketsCount(0), isCached(false)
  {
  }

  ElementIterator Begin(IndexType bucketIndex) const
  {
    return rangeBegins[bucketIndex] < elementIndices.size() ? &elementIndices[rangeBegins[bucketIndex]] : 0;
  }

  ElementIterator End(IndexType bucketIndex) const
  {
    return Begin(bucketIndex) + rangeEnds[bucketIndex] - rangeBegins[bucketIndex];
  }

  void Initialize(IndexType bucketsCount, IndexType storageSize)
  {
    this->bucketsCount = bucketsCount;
    assert(bucketsCount > 0);

    elementLocations.resize(storageSize);
    elementIndices.resize(storageSize);
    bucketIndices.resize(storageSize);
    rangeBegins.resize(bucketsCount);
    rangeEnds.resize(bucketsCount);


    std::fill(rangeBegins.begin(), rangeBegins.end(), 0);
    std::fill(rangeEnds.begin(),   rangeEnds.end(),   0);

    #pragma omp parallel for
    for (int elementIndex = 0; elementIndex < int(storageSize); ++elementIndex)
    {
      elementIndices[elementIndex] = elementIndex;
      elementLocations[elementIndices[elementIndex]] = elementIndex;
    }

    rangeBegins.back() = 0;
    rangeEnds.back() = storageSize;
    isCached = false;
  }

  void Cache()
  {
    isCached = true;
    #pragma omp parallel for
    for (int bucketIndex = 0; bucketIndex < int(bucketsCount); ++bucketIndex)
    {
      for (IndexType elementLocation = rangeBegins[bucketIndex]; elementLocation < rangeEnds[bucketIndex]; ++elementLocation)
      {
        bucketIndices[elementIndices[elementLocation]] = bucketIndex;
      }
    }
  }

  IndexType operator[](IndexType elementIndex) const
  {
    return GetBucket(elementIndex);
  }

  void UpdateBucket(IndexType elementIndex, IndexType newBucket)
  {
    IndexType currentBucket = GetBucket(elementIndex);
    while (currentBucket < newBucket)
    {
      MoveToNextBucket(elementIndex);
      ++currentBucket;
    }

    while (currentBucket > newBucket)
    {
      MoveToPreviousBucket(elementIndex);
      --currentBucket;
    }
  }

  std::vector<IndexType> elementIndices;

private:
  IndexType bucketsCount;
  std::vector<IndexType> elementLocations;
  std::vector<IndexType> bucketIndices;
  std::vector<IndexType> rangeBegins;
  std::vector<IndexType> rangeEnds;
  bool isCached;

  IndexType GetBucket(IndexType elementIndex) const
  {
    if (isCached)
    {
      return bucketIndices[elementIndex];
    }

    IndexType bucketIndex = 0;
    while (bucketIndex < bucketsCount)
    {
      if (rangeBegins[bucketIndex] <= elementLocations[elementIndex] && elementLocations[elementIndex] < rangeEnds[bucketIndex])
      {
        break;
      }
      ++bucketIndex;
    }
    assert(bucketIndex < bucketsCount);
    return bucketIndex;
  }

  void MoveToNextBucket(IndexType elementIndex)
  {
    IndexType currentBucket = GetBucket(elementIndex);
    assert(currentBucket + 1 < bucketsCount);

    IndexType currentLocation = elementLocations[elementIndex];
    Swap(currentLocation, rangeEnds[currentBucket] - 1);

    --rangeBegins[currentBucket + 1];
    --rangeEnds[currentBucket];
  }

  void MoveToPreviousBucket(IndexType elementIndex)
  {
    IndexType currentBucket = GetBucket(elementIndex);
    assert(currentBucket > 0);

    IndexType currentLocation = elementLocations[elementIndex];
    Swap(currentLocation, rangeBegins[currentBucket]);

    ++rangeBegins[currentBucket];
    ++rangeEnds[currentBucket - 1];
  }

  void Swap(IndexType firstLocation, IndexType secondLocation)
  {
    std::swap(elementIndices[firstLocation], elementIndices[secondLocation]);
    elementLocations[elementIndices[firstLocation]] = firstLocation;
    elementLocations[elementIndices[secondLocation]] = secondLocation;
  }
};

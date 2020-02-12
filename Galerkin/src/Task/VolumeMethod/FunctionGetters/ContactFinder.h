#pragma once

#include "../../../Maths/Collisions.h"
#include "../../GeomMesh/AABBTree.h"

template <typename ContactFinder>
struct ContactProcessor
{
  typedef typename ContactFinder::MeshTypeT MeshType;
  typedef typename MeshType::Space          Space;
  typedef typename MeshType::SystemT        System;
  typedef typename System::ValueType        ValueType;
  SPACE_TYPEDEFS

  const static int dimsCount = MeshType::dimsCount;

  ContactProcessor(ContactFinder* contactFinder,
    IndexType cellIndex, IndexType* contactedCells, IndexType* contactedCellsCount, IndexType maxContactCellsCount):
    contactFinder(contactFinder), cellIndex(cellIndex), contactedCells(contactedCells), contactedCellsCount(contactedCellsCount), maxContactCellsCount(maxContactCellsCount)
  {
    (*contactedCellsCount) = 0;
  }

  void operator()(int nodeIndex)
  {
    IndexType neighbourCellIndex = contactFinder->mesh->aabbTree.GetUserData(nodeIndex);
    if (contactFinder->TryCell(cellIndex, neighbourCellIndex))
    {
      assert(!(contactedCellsCount && *contactedCellsCount + 1 >= maxContactCellsCount)); //max number of contacts per cell exceeded, usually means something went terribly wrong
      if (contactedCellsCount && *contactedCellsCount + 1 < maxContactCellsCount)
      {
        if (contactedCells)
          contactedCells[*contactedCellsCount] = neighbourCellIndex;

        ++(*contactedCellsCount);
      }
    }
  }

  ContactFinder* contactFinder;
  IndexType cellIndex;
  IndexType* contactedCells; // pointer to array of colliding cells
  IndexType* contactedCellsCount;
  IndexType maxContactCellsCount;
};


template<typename MeshType>
struct ContactFinder
{
  typedef typename MeshType::Space          Space;
  typedef typename MeshType::SystemT        System;
  typedef MeshType MeshTypeT;
  typedef typename System::ValueType        ValueType;
  typedef typename System::MediumParameters MediumParameters;
  SPACE_TYPEDEFS

  typedef ContactProcessor< ContactFinder<MeshType> > ContactProcessorT;

  const static int dimsCount = MeshType::dimsCount;

  ContactFinder(
    MeshType* mesh,
    IndexType cellIndex):
    mesh(mesh),
    cellIndex(cellIndex)
  {
    mesh->GetFixedCellIndices(cellIndex, cellIndices);
    mesh->GetCellVertices(cellIndices, cellVertices);
  }

  MeshType* mesh;
  IndexType cellIndex;
  IndexType cellIndices[Space::NodesPerCell];
  Vector cellVertices[Space::NodesPerCell];

  void Find(IndexType* collidedCells, IndexType* collidedCellsCount, IndexType maxCollidedCellsCount)
  {
    AABB aabb = mesh->GetCellAABB(cellIndex);
    ContactProcessorT contactProcessor(this, cellIndex, collidedCells, collidedCellsCount, maxCollidedCellsCount);
    mesh->aabbTree.template FindCollisions<ContactProcessorT>(aabb, contactProcessor);
  }

  bool TryCell(IndexType cellIndex, IndexType neighbourCellIndex)
  {
    if (cellIndex == neighbourCellIndex) return false;

    Vector neighbourCellVertices[Space::NodesPerCell];
    mesh->GetCellVertices(neighbourCellIndex, neighbourCellVertices);

    return CellsCollide(cellVertices, neighbourCellVertices, mesh->collisionWidth);
    // return CellsCollide(cellVertices, neighbourCellVertices); // for 2d only
  }
};

template <typename TreeContactFinder>
struct TreeContactProcessor
{
  typedef typename TreeContactFinder::MeshTypeT MeshType;
  typedef typename MeshType::Space          Space;
  typedef typename MeshType::SystemT        System;
  typedef typename System::ValueType        ValueType;
  SPACE_TYPEDEFS

  const static int dimsCount = MeshType::dimsCount;

  TreeContactProcessor(TreeContactFinder* contactFinder) : contactFinder(contactFinder)
  {}

  void operator()(int nodeIndex0, int nodeIndex1)
  {
    const IndexType cellIndex0 = contactFinder->mesh->aabbTree.GetUserData(nodeIndex0);
    const IndexType cellIndex1 = contactFinder->mesh->aabbTree.GetUserData(nodeIndex1);

    if (contactFinder->TryCell(cellIndex0, cellIndex1))
    {
      IndexType& cellsCount0 = contactFinder->mesh->cellsCollided[cellIndex0].count;
      contactFinder->mesh->cellsCollided[cellIndex0].cells[cellsCount0++] = cellIndex1;

      IndexType& cellsCount1 = contactFinder->mesh->cellsCollided[cellIndex1].count;
      contactFinder->mesh->cellsCollided[cellIndex1].cells[cellsCount1++] = cellIndex0;
    }
  }

  TreeContactFinder* contactFinder;
};

template<typename MeshType>
struct TreeContactFinder
{
  typedef typename MeshType::Space          Space;
  typedef typename MeshType::SystemT        System;
  typedef MeshType MeshTypeT;
  typedef typename System::ValueType        ValueType;
  typedef typename System::MediumParameters MediumParameters;

  SPACE_TYPEDEFS

  typedef AABBTree<Space, IndexType> AABBTreeT;
  typedef TreeContactProcessor< TreeContactFinder<MeshType> > TreeContactProcessorT;

  const static int dimsCount = MeshType::dimsCount;

  TreeContactFinder(
    MeshType* mesh) :
    mesh(mesh)
  {
    for (IndexType cellIndex = 0; cellIndex < mesh->cells.size(); ++cellIndex)
    {
      mesh->cellsCollided[cellIndex].count = 0;
    }
  }

  MeshType* mesh;

  void Find(const AABBTreeT& other)
  {
    TreeContactProcessorT contactProcessor(this);
    mesh->aabbTree.template FindTreeCollisions< AABBTreeT, TreeContactProcessorT >(other, contactProcessor);
  }

  bool TryCell(IndexType cellIndex0, IndexType cellIndex1)
  {
    if (cellIndex0 == cellIndex1) return false;

    Vector cellVertices0[Space::NodesPerCell];
    mesh->GetCellVertices(cellIndex0, cellVertices0);

    Vector cellVertices1[Space::NodesPerCell];
    mesh->GetCellVertices(cellIndex1, cellVertices1);

    return CellsCollide(cellVertices0, cellVertices1, mesh->collisionWidth);
  }
};



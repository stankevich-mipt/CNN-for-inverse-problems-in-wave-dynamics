#pragma once

template <typename GhostCellFunctionGetter>
struct CollisionProcessor
{
  typedef typename GhostCellFunctionGetter::MeshTypeT MeshType;
  typedef typename MeshType::Space   Space;

  typedef typename MeshType::SystemT     System;
  typedef typename System::ValueType     ValueType;
  SPACE_TYPEDEFS

  const static int dimsCount = MeshType::dimsCount;

  CollisionProcessor(GhostCellFunctionGetter* functionGetter,
    const Vector& testPoint, Scalar* values):
    functionGetter(functionGetter), testPoint(testPoint), values(values), found(false), collidedCellIndex(IndexType(-1))
  {}
  void operator()(int nodeIndex)
  {
    if (!found)
    {
      IndexType neighbourCellIndex = functionGetter->mesh->aabbTree.GetUserData(nodeIndex);
      if (functionGetter->TryGhostCell(testPoint, neighbourCellIndex, values))
      {
        found = true;
        collidedCellIndex = neighbourCellIndex;
      }
    }
  }

  GhostCellFunctionGetter* functionGetter;
  Vector testPoint;
  Scalar* values;
  bool found;
  IndexType collidedCellIndex;
};

template<typename MeshType>
struct GhostCellFunctionGetter
{
  typedef typename MeshType::Space   Space;

  typedef typename MeshType::SystemT     System;
  typedef MeshType MeshTypeT;
  typedef typename System::ValueType     ValueType;
  typedef typename System::MediumParameters MediumParameters;
  SPACE_TYPEDEFS

  typedef CollisionProcessor< GhostCellFunctionGetter<MeshType> > CollisionProcessorT;

  const static int dimsCount = MeshType::dimsCount;

  enum GetterType
  {
    Solution, MediumParams
  };

  GhostCellFunctionGetter(
    MeshType* mesh,
    IndexType cellIndex,
    Scalar currTime,
    GetterType getterType):
      cellIndex(cellIndex),
      mesh(mesh),
      currTime(currTime),
      getterType(getterType)
  {
  }

  IndexType cellIndex;
  MeshType* mesh;
  Scalar currTime;
  GetterType getterType;

  IndexType operator()(const Vector& testPoint, Scalar* values)
  {
    switch (getterType)
    {
      case Solution: std::fill(values, values + dimsCount, 0); break;
      case MediumParams: std::fill(values, values + MediumParameters::ParamsCount, 0); break; // lambda, mju, invRho 
    }

    CollisionProcessorT collisionProcessor(this, testPoint, values);
    mesh->aabbTree.template FindCollisions<CollisionProcessorT>(AABB(testPoint, testPoint), collisionProcessor);

    return collisionProcessor.collidedCellIndex;
  }

  IndexType operator()(const Vector& testPoint, IndexType* const contactCells, IndexType contactCellsCount, Scalar* values)
  {
    switch (getterType)
    {
      case Solution: std::fill(values, values + dimsCount, 0); break;
      case MediumParams: std::fill(values, values + MediumParameters::ParamsCount, 0); break; // lambda, mju, invRho 
    }

    IndexType collidedCellIndex = IndexType(-1);
    for (IndexType cellNumber = 0; cellNumber < contactCellsCount; ++cellNumber)
    {
      IndexType neighbourCellIndex = contactCells[cellNumber];
      if (TryGhostCell(testPoint, neighbourCellIndex, values))
      {
        collidedCellIndex = neighbourCellIndex;
        break;
      }
    }
    return collidedCellIndex;
  }

  bool TryGhostCell(const Vector& testPoint, int neighbourCellIndex,
    Scalar* values)
  {
    if (cellIndex == neighbourCellIndex) return false;

    IndexType neighbourCellIndices[Space::NodesPerCell];
    mesh->GetCellIndices(neighbourCellIndex, neighbourCellIndices);

    Vector neighbourCellVertices[Space::NodesPerCell];
    mesh->GetCellVertices(neighbourCellIndices, neighbourCellVertices);

    //point is inside a colliding cell, counting it as a contact
    if (PointInCell<Scalar>(neighbourCellVertices, testPoint))
    {
      switch (getterType)
      {
        case Solution: 
        {
          ValueType neighbourSolution = mesh->GetCellSolution(neighbourCellIndex, testPoint);
          for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
          {
            values[valueIndex] = neighbourSolution.values[valueIndex];
          }
        }
        break;
        case MediumParams:
          for (IndexType paramIndex = 0; paramIndex < MediumParameters::ParamsCount; ++paramIndex)
          {
            values[paramIndex] = mesh->cellMediumParameters[neighbourCellIndex].params[paramIndex];
          }
        break;
      }
      return true;
    }
    return false;
  }
};

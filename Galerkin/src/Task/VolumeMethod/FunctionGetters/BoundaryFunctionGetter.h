#pragma once
#include "../../../Maths/Spaces.h"
#include "../../ElasticSystem/BoundaryInfos.h"

template <typename MeshType>
class BoundaryFunctionGetter
{
public:
  typedef typename MeshType::Space Space;
  SPACE_TYPEDEFS

  typedef typename MeshType::SystemT System;
  typedef typename System::ValueType ValueType;

  BoundaryFunctionGetter(BoundaryInfoFunctor<Space>* functor,
    const Scalar time,
    const MeshType* const volumeMesh,
    IndexType cellIndex, IndexType boundaryFaceNumber):
    functor(functor),
    time(time),
    volumeMesh(volumeMesh),
    cellIndex(cellIndex)
  {
    volumeMesh->GetGhostCellVertices(cellIndex, boundaryFaceNumber, ghostPoints);
    BuildExternalNormal(cellIndex, boundaryFaceNumber);
  }

  void operator()(const Vector& refPoint, Scalar* values)
  {
    if (functor)
    {
      ValueType valueType = volumeMesh->GetRefCellSolution(cellIndex, refPoint);
      functor->SetCurrentVelocity(valueType.GetVelocity());
      Vector globalPoint = volumeMesh->RefToGlobalVolumeCoords(refPoint, ghostPoints);
      functor->operator()(globalPoint, externalNormal, time, values);
    } else
    {
      std::fill_n(values, MeshType::dimsCount, Scalar(0.0));
    }
  }

  void BuildExternalNormal(IndexType cellIndex, IndexType boundaryFaceNumber)
  {
    BuildExternalNormal(Overload<Space>(), cellIndex, boundaryFaceNumber);
  }

private:
  void BuildExternalNormal(Overload<Space2>, IndexType cellIndex, IndexType boundaryEdgeNumber)
  {
    externalNormal = volumeMesh->GetEdgeExternalNormal(cellIndex, boundaryEdgeNumber);
  }

  void BuildExternalNormal(Overload<Space3>, IndexType cellIndex, IndexType boundaryFaceNumber)
  {
    externalNormal = volumeMesh->GetFaceExternalNormal(cellIndex, boundaryFaceNumber);
  }

  BoundaryInfoFunctor<Space>* functor;
  Scalar time;
  const MeshType* const volumeMesh;
  Vector ghostPoints[Space::NodesPerCell];
  Vector externalNormal;
  IndexType cellIndex;
};

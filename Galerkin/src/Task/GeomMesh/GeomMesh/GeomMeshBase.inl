GeomMeshBase<Space2>::EdgeLocation::EdgeLocation(): cellIndex(-1), edgeNumber(-1)
{
}

GeomMeshBase<Space2>::EdgeLocation::EdgeLocation(IndexType cellIndex, IndexType edgeNumber):
  cellIndex(cellIndex), edgeNumber(edgeNumber)
{
}

bool GeomMeshBase<Space2>::EdgeLocation::IsNull() const
{
  return cellIndex == IndexType(-1) || edgeNumber == IndexType(-1);
}

bool GeomMeshBase<Space2>::EdgeLocation::operator == (const EdgeLocation& other) const
{
  return cellIndex == other.cellIndex && edgeNumber == other.edgeNumber;
}

bool GeomMeshBase<Space2>::EdgeLocation::operator != (const EdgeLocation& other) const
{
  return !this->operator==(other);
}

GeomMeshBase<Space2>::EdgeLocationPair::EdgeLocationPair()
{
}

GeomMeshBase<Space2>::EdgeLocationPair::EdgeLocationPair(
  const EdgeLocation& firstEdge, const EdgeLocation& secondEdge)
{
  edges[0] = firstEdge;
  edges[1] = secondEdge;
}

GeomMeshBase<Space3>::FaceLocation::FaceLocation(): cellIndex(-1), faceNumber(-1)
{
}

GeomMeshBase<Space3>::FaceLocation::FaceLocation(IndexType cellIndex, IndexType faceNumber):
  cellIndex(cellIndex), faceNumber(faceNumber)
{
}

bool GeomMeshBase<Space3>::FaceLocation::IsNull() const
{
  return cellIndex == IndexType(-1) || faceNumber == IndexType(-1);
}

bool GeomMeshBase<Space3>::FaceLocation::operator == (const FaceLocation& other) const
{
  return cellIndex == other.cellIndex && faceNumber == other.faceNumber;
}

bool GeomMeshBase<Space3>::FaceLocation::operator != (const FaceLocation& other) const
{
  return !this->operator==(other);
}

GeomMeshBase<Space3>::FaceLocationPair::FaceLocationPair()
{
}

GeomMeshBase<Space3>::FaceLocationPair::FaceLocationPair(
  const FaceLocation& firstFace, const FaceLocation& secondFace)
{
  faces[0] = firstFace;
  faces[1] = secondFace;
}


#pragma once

#include "MeshIOBase.h"
#include "../../../../IO/Vtk/MeshVtkWriter.h"
#include "../../../../Utils/Utils.h"

template <typename Space>
struct MeshIO: public MeshIOBase<Space>
{
  SPACE_TYPEDEFS
  MeshIO();
  virtual ~MeshIO() = default;
  virtual void Save(const std::string& fileName, IO::FileType fileType = IO::Binary);
  virtual void Load(const std::string& fileName, IO::FileType fileType = IO::Binary);

  // file should be opened
  void Save(std::fstream& file, IO::FileType fileType);
  void Load(std::fstream& file, IO::FileType fileType);

  std::vector<Vector>    vertices;
  std::vector<IndexType> indices;

  IndexType contactTypesCount;
  IndexType boundaryTypesCount;

  std::vector<Vector> detectorsPositions;

  IndexType GetCellsCount() const;
  Vector GetCellCenter(IndexType cellIndex) const;
  
  bool operator==(const MeshIO<Space>& other) const;

  void SaveContacts(const std::string& vtkFileName);
  void SaveBoundaries(const std::string& vtkFileName);
  void SaveDetectors(const std::string& vtkFileName);

private:
  void SaveContacts(std::fstream& file, IO::FileType fileType);
  void SaveBoundaries(std::fstream& file, IO::FileType fileType);
  void LoadContacts(std::fstream& file, IO::FileType fileType);
  void LoadBoundaries(std::fstream& file, IO::FileType fileType);
};

#include "MeshIO.inl"
#include "MeshIO2.inl"
#include "MeshIO3.inl"

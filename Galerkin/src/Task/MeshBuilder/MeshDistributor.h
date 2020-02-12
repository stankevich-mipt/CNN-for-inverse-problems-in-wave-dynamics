#pragma once

#include <vector>
#include "../GeomMesh/MeshIO/Local/MeshIO.h"
#include "../../Maths/Spaces.h"
#include <metisbin.h>
#include <numeric>

template <typename Space>
class MeshDestributor
{
public:
  SPACE_TYPEDEFS
  MeshDestributor()
  {
    printf("Distributing cells between domains\n");
  }

  virtual ~MeshDestributor() = default;

  virtual void Distribute(const MeshIO<Space>* const mesh, IndexType domainsCount, 
    const std::vector<IndexType>& cellsComputationalCosts, const std::string& algorithmName, std::vector<IndexType>* const cellsColors) const = 0;
};

template <typename Space>
class MetisDistributor: public MeshDestributor<Space>
{
public:
  SPACE_TYPEDEFS

  void Distribute(const MeshIO<Space>* const mesh, IndexType domainsCount, 
    const std::vector<IndexType>& cellsComputationalCosts, 
    const std::string& algorithmName, std::vector<IndexType>* const cellsColors) const override
  {
    cellsColors->resize(mesh->GetCellsCount(), 0);
    if (domainsCount == 1) return;

    int ne = mesh->GetCellsCount();
    int nn = mesh->vertices.size();
    idx_t* eptr =static_cast<idx_t*>(malloc((ne + 1) * sizeof(idx_t)));
    for (IndexType cellIndex = 0; cellIndex < IndexType(ne + 1); ++cellIndex)
    {
      eptr[cellIndex] = Space::NodesPerCell * cellIndex;
    }
    idx_t* eind = static_cast<idx_t*>(malloc(sizeof(idx_t) * mesh->indices.size()));
    std::copy(mesh->indices.begin(), mesh->indices.end(), eind);

    /* "ncommon" specifies the number of common nodes that two elements must have in 
       order to put an edge between them in the dual graph */
    idx_t ncommon = Space::Dimension;
    idx_t nparts = domainsCount;
    idx_t objval = 0;

    idx_t* epart = static_cast<idx_t*>(malloc(ne * sizeof(idx_t)));
    idx_t* npart = static_cast<idx_t*>(malloc(nn * sizeof(idx_t)));
    idx_t* vwgt  = static_cast<idx_t*>(malloc(ne * sizeof(idx_t)));
    for (IndexType cellIndex = 0; cellIndex < IndexType(ne); ++cellIndex)
    {
      vwgt[cellIndex] = static_cast<idx_t>(cellsComputationalCosts[cellIndex]);
    }
    
    int status = 0;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
    options[METIS_OPTION_CONTIG] = 1; // Forces contiguous partitions

    // METIS_PartMeshDual is better
    if (algorithmName == "METIS_PartMeshDual")
      status = ::METIS_PartMeshDual(&ne, &nn, eptr, eind, vwgt, nullptr, &ncommon, &nparts, nullptr, options, &objval, epart, npart);
    else if (algorithmName == "METIS_PartMeshNodal")
      status = ::METIS_PartMeshNodal(&ne, &nn, eptr, eind, vwgt, nullptr, &nparts, nullptr, options, &objval, epart, npart);
    else 
    {
      std::cerr << "Unknown partitioning algorithm\n";
      assert(0);
    }

    if (status == METIS_OK)
    {
      std::copy(epart, epart + ne, cellsColors->begin());
    } else {
      std::cerr << "Metis can`t distribute mesh\n";
    }

    std::cout << "Domains count = " << domainsCount 
              << "; Total communication value = " << objval << std::endl;

    free(epart);
    free(npart);
    free(eind);
    free(eptr);
    free(vwgt);
  }
};

template <typename Space>
class SimpleRectDistributor;

template <>
class SimpleRectDistributor<Space2>: public MeshDestributor<Space2>
{
public:
  using Space = Space2;
  SPACE_TYPEDEFS

  void Distribute(const MeshIO<Space>* const mesh, IndexType domainsCount, const std::vector<IndexType>&,
    const std::string&, std::vector<IndexType>* const cellsColors) const override
  {
    AABB boundingBox;
    for (IndexType nodeIndex = 0; nodeIndex < mesh->vertices.size(); ++nodeIndex)
    {
      boundingBox.Expand(mesh->vertices[nodeIndex]);
    }

    std::vector<Scalar> horizonatalStrips;
    std::vector<Scalar> verticalStrips;

    IndexVector proportions(1, domainsCount);

    for (IndexType horizontalStripsCount = 2; horizontalStripsCount <= domainsCount; ++horizontalStripsCount)
    {
      if (domainsCount % horizontalStripsCount == 0)
      {
        IndexType verticalStripsCount = domainsCount / horizontalStripsCount;
        Vector boxSize(boundingBox.boxPoint2.x - boundingBox.boxPoint1.x, 
                       boundingBox.boxPoint2.y - boundingBox.boxPoint1.y);
        if (fabs(boxSize.y * Scalar(proportions.x) - boxSize.x * Scalar(proportions.y)) > 
            fabs(boxSize.y * Scalar(horizontalStripsCount) - boxSize.x * Scalar(verticalStripsCount)))
        {
          proportions = IndexVector(horizontalStripsCount, verticalStripsCount);
        }
      }
    }

    horizonatalStrips.push_back(-std::numeric_limits<Scalar>::max() * Scalar(0.5));
    Scalar horizonatalStep = (boundingBox.boxPoint2.y - boundingBox.boxPoint1.y) / proportions.y;
    for (IndexType stripIndex = 1; stripIndex < proportions.y; ++stripIndex)
    {
      horizonatalStrips.push_back(boundingBox.boxPoint1.y + horizonatalStep * Scalar(stripIndex));
    }
    horizonatalStrips.push_back(std::numeric_limits<Scalar>::max() * Scalar(0.5));

    verticalStrips.push_back(-std::numeric_limits<Scalar>::max() * Scalar(0.5));
    Scalar verticalStep = (boundingBox.boxPoint2.x - boundingBox.boxPoint1.x) / proportions.x;
    for (IndexType stripIndex = 1; stripIndex < proportions.x; ++stripIndex)
    {
      verticalStrips.push_back(boundingBox.boxPoint1.x + verticalStep * Scalar(stripIndex));
    }
    verticalStrips.push_back(std::numeric_limits<Scalar>::max() * Scalar(0.5));

    cellsColors->resize(mesh->indices.size() / Space2::NodesPerCell);
    for (IndexType cellIndex = 0; cellIndex < cellsColors->size(); ++cellIndex)
    {
      Vector2 centre(0, 0);
      for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerCell; ++nodeNumber)
      {
        IndexType nodeIndex = mesh->indices[Space2::NodesPerCell * cellIndex + nodeNumber];
        centre += mesh->vertices[nodeIndex];
      }
      centre /= Scalar(Space2::NodesPerCell);

      IndexVector domainPosition = IndexVector::zero();

      while (centre.x >    verticalStrips[domainPosition.x + 1]) domainPosition.x++;
      while (centre.y > horizonatalStrips[domainPosition.y + 1]) domainPosition.y++;

      (*cellsColors)[cellIndex] = domainPosition.x + domainPosition.y * proportions.x;
    }
    printf("\n\n");
  }
};

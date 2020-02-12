#pragma once

//#include "glib/LinearSystemSolvers.h"
#include "ContactCorrectors.h"


template <typename SpaceT>
class ContactGcCalculator3D
{
  typedef typename SpaceT::Tensor3 Tensor3;
  typedef typename SpaceT::ElasticProperties ElasticProperties;
  typedef typename SpaceT::Elastic Elastic;
  typedef typename SpaceT::Vector3 Vector3;
  typedef typename SpaceT::Scalar Real;

  typedef typename SpaceT::IndexType  IndexType;
public:

  ContactGcCalculator3D(){ }

  bool randomizeSplitAxes() const { return d_randomizeSplitAxes; }
  void randomizeSplitAxes(bool i_randomizeSplitAxes)
    {d_randomizeSplitAxes = i_randomizeSplitAxes; }

  Elastic getMaterial(const size_t i_materialId) const { return d_cachedMaterials[i_materialId].d_material; }
  void setMaterial(const ElasticProperties & i_material, size_t i_materialId, size_t stepMultiplier);

public:
  void generateSplitVectors();
  void initiateSplitStep(size_t i_splitPhase, Real i_step, size_t i_globalStepNumber);

  void SetBoundarySettings(const BoundarySettings boundarySettings, const IndexType typeIndex);
  void SetContactSettings (const ContactSettings contactSettings, const IndexType typeIndex);

  template <typename It>
  void transformWaves(It i_begin, It i_end) const;

  template <typename It>
  void transformDirect(It i_begin, It i_end) const;

  // compute next nodes' values;
  // it is assumed that waves array is stored in nodes (instead of elastic)
  template <typename It, typename Reconstructor>
  void computeInnerNodesWaves(It i_begin, It i_end, const Reconstructor & i_r) const;
  template <typename It, typename Reconstructor>
  void computeOuterNodesWaves(It i_begin, It i_end, const Reconstructor & i_r) const;

  template <typename It, typename Reconstructor>
  void computeInnerNodesDirect(It i_begin, It i_end, const Reconstructor & i_r) const;
  template <typename It, typename Reconstructor>
  void computeOuterNodesDirect(It i_begin, It i_end, const Reconstructor & i_r) const;

  Vector3 getSplitDirection();

  template <typename It>
  void copyNextToCurr(It i_begin, It i_end) const;

private:
  enum Wave3D
  {
    WAVE_3D_L1, WAVE_3D_L2, WAVE_3D_L3,
    WAVE_3D_R3, WAVE_3D_R2, WAVE_3D_R1,
    WAVE_3D_COUNT
  };

  // multiplies given elastic (v, s) on six eigenrows (Omega matrix)
  // and places the result in waves (elastic and waves array overlap in memory)
  void computeWaves(const Elastic & i_e, const size_t i_materialId, Real o_waves[WAVE_3D_COUNT]) const;
  // delta1 = w1^n(-c1*tau) - w1^n(0)
  // delta2 = w3^n(-c2*tau) - w3^n(0)
  // delta3 = w5^n(-c2*tau) - w5^n(0)
  Elastic incLeftWaves (Real i_delta1, Real i_delta2, Real i_delta3, const size_t i_materialId) const;
  // delta1 = w2^n(+c1*tau) - w2^n(0)
  // delta2 = w4^n(+c2*tau) - w4^n(0)
  // delta3 = w6^n(+c2*tau) - w6^n(0)
  Elastic incRightWaves(Real i_delta1, Real i_delta2, Real i_delta3, const size_t i_materialId) const;

  Elastic incLeftDirect(const Elastic &lc1, const Elastic &lc2, const size_t i_materialId) const;
  Elastic incRightDirect(const Elastic &rc1, const Elastic &rc2, const size_t i_materialId) const;
  Elastic incDirect(const Elastic &lc1, const Elastic &lc2, const Elastic &rc1, const Elastic &rc2, const size_t i_materialId) const;

  void correctFreeBoundary(Elastic &point, const Vector3 &normal, ElasticProperties material, bool isLeftSide) const;


private:
  bool d_randomizeSplitAxes;

  Vector3 d_n;        // split direction of the current split phase
  Vector3 d_n1, d_n2; // (n, n1, n2) form orthogonal coordinate frame

  Tensor3 d_N;       // d_n times d_n
  Tensor3 d_N1;      // 1/2 (tensor-product(n, n1) + tensor-product(n1, n))
  Tensor3 d_N2;      // 1/2 (tensor-product(n, n2) + tensor-product(n2, n))

  Vector3 basisVectors[3];

  size_t globalStepNumber;

  std::vector<ContactSettings> contactSettings;
  std::vector<BoundarySettings> boundarySettings;

  struct CachedMaterial
  {
    // cached values to speed up calculations
    Real d_rrhoc1;     // 1 / (rho * c1)
    Real d_rrhoc2;     // 1 / (rho * c2)
    Real d_rhoc2;      // rho * c2
    Real d_rhoc3;      // rho * c3
    Real d_rhoc1mc3;   // rho * (c1 - c3)
    Vector3 d_probeC1;  // c1 * step * n
    Vector3 d_probeC2;  // c2 * step * n
    size_t stepMultiplier;
    ElasticProperties d_material;
  };

  std::vector<CachedMaterial> d_cachedMaterials;
};


// explicit instantiation
//template class ContactGcCalculator3D<float>;
//template class ContactGcCalculator3D<double>;


#include "ContactGcCalculator3D.inl"

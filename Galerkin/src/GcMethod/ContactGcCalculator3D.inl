//#include "gcmethod/GcCalculator3D.h"



template <typename Space>
void ContactGcCalculator3D<Space>::setMaterial(const ElasticProperties & i_material, const size_t i_materialId, const size_t stepMultiplier)
{
  if(i_materialId + 1 > d_cachedMaterials.size()) d_cachedMaterials.resize(i_materialId + 1);

  d_cachedMaterials[i_materialId].d_material = i_material;
  Real rho = Real(1.0) / i_material.invRho;
  Real c1 = i_material.c1();
  Real c2 = i_material.c2();
  Real c3 = i_material.c3();
  d_cachedMaterials[i_materialId].d_rrhoc1 = 1 / (rho * c1);
  d_cachedMaterials[i_materialId].d_rrhoc2 = 1 / (rho * c2);
  d_cachedMaterials[i_materialId].d_rhoc2 = rho * c2;
  d_cachedMaterials[i_materialId].d_rhoc3 = rho * c3;
  d_cachedMaterials[i_materialId].d_rhoc1mc3 = rho * (c1 - c3);
  d_cachedMaterials[i_materialId].stepMultiplier = stepMultiplier;
}

template <typename Space>
template <typename It>
void ContactGcCalculator3D<Space>::transformWaves(It i_begin, It i_end) const
{
  const long int n = (long int)(i_end - i_begin);
  //#pragma omp parallel for
  for (long int i = 0; i < n; ++i)
  {
    (i_begin + i).nextValue() = (i_begin + i).currValue(); //DEBUG
    (i_begin + i).tempValue() = (i_begin + i).currValue();

    Elastic & e = (i_begin + i).currValue();
    computeWaves(e, (i_begin + i).GetSubmeshIndex(), &e(0));

    /*Elastic & nextElastic = (i_begin + i).nextValue();
    Tensor3 N11 = SymmetricTensorProduct(d_n1, d_n1) * Real(0.5);
    Tensor3 N12 = SymmetricTensorProduct(d_n1, d_n2) * Real(0.5);
    Tensor3 N22 = SymmetricTensorProduct(d_n2, d_n2) * Real(0.5);

    Real lambda = d_cachedMaterials[(i_begin + i).GetSubmeshIndex()].d_material.lambda;
    Real mu = d_cachedMaterials[(i_begin + i).GetSubmeshIndex()].d_material.mju;

    Tensor3 I(Real(1.0));

    Real w7 = DoubleConvolution(nextElastic.GetTension(), N12);
    Real w8 = DoubleConvolution(nextElastic.GetTension(), N11 - N22);
    Real w9 = DoubleConvolution(nextElastic.GetTension(), (N11 + N22 - d_N * (Real(2.0) * lambda * (lambda + Real(2.0) * mu))));


    nextElastic.s += Real(2.0) * w7 * N12;
    nextElastic.s += Real(0.5) * w8 * (N11 - N22);
    nextElastic.s += Real(0.5) * w9 * (I - d_N);*/
  }
}

template <typename Space>
template <typename It>
void ContactGcCalculator3D<Space>::transformDirect(It i_begin, It i_end) const
{
  const long int n = (long int)(i_end - i_begin);
  #pragma omp parallel for
  for (long int i = 0; i < n; ++i)
  {
    (i_begin + i).nextValue() = (i_begin + i).currValue();
    Elastic & e = (i_begin + i).currValue();
  }
}

template <typename Space>
void ContactGcCalculator3D<Space>::correctFreeBoundary(Elastic &point, const Vector3 &normal, ElasticProperties material, bool isLeftSide) const
{
  Vector3 z = point.GetTension()(normal);
  Real wmega1 = Real(2.0) * (d_n * normal) * (d_n * z) - (normal * z);

  Real divider = (material.c1() + material.c3()) * sqr((d_n * normal)) - material.c3();
  if(fabs(divider) > std::numeric_limits<Real>::epsilon())
  {
    wmega1 /= divider;
  }else
  {
    return;
  }

  Vector3 b;
  if(fabs((d_n * normal)) > std::numeric_limits<Real>::epsilon())
  {
    b = (z - wmega1 * material.c3() * normal) / (material.c2() * (d_n * normal));
  }else
  {
    return;
  }
  Real sideMultiplier = !isLeftSide ? Real(1.0) : Real(-1.0);

  point.SetVelocity(point.GetVelocity() - sideMultiplier * material.invRho * ((wmega1 - b * d_n) * d_n + b));
  point.SetTension (point.GetTension () - d_N * ((material.c1() - material.c3()) * wmega1 - Real(2.0) * material.c2() * (d_n * b)));
  point.SetTension (point.GetTension () - Tensor3(material.c3() * wmega1));
  point.SetTension (point.GetTension () - Tensor3(material.c2() * (d_n * b)));
}






template <typename Space>
template <typename It, typename Reconstructor>
void ContactGcCalculator3D<Space>::computeInnerNodesWaves(It i_begin, It i_end,
  const Reconstructor & i_r) const
{
  const long int n = (long int)(i_end - i_begin);
  //#pragma omp parallel for
  for (long int i = 0; i < n; i++)
  {
    It it = i_begin + i;
    const Elastic & currElastic = it.currValue();
    Elastic & nextElastic = it.nextValue();
    Vector3 pos = it.currPos();



    size_t negCell, posCell;
    it.findCellAlongDir(d_n, negCell, posCell);

    if ((negCell != (size_t)(-1)) && (posCell != (size_t)(-1)))
    {
      {
        const Real DELTA_1 =
          i_r.template interpolateInCell<WAVE_3D_L1>(negCell, pos - d_cachedMaterials[it.GetSubmeshIndex()].d_probeC1) -
          currElastic(WAVE_3D_L1);
        const Real DELTA_2 =
          i_r.template interpolateInCell<WAVE_3D_L2>(negCell, pos - d_cachedMaterials[it.GetSubmeshIndex()].d_probeC2) -
          currElastic(WAVE_3D_L2);
        const Real DELTA_3 =
          i_r.template interpolateInCell<WAVE_3D_L3>(negCell, pos - d_cachedMaterials[it.GetSubmeshIndex()].d_probeC2) -
          currElastic(WAVE_3D_L3);
        nextElastic += incLeftWaves(DELTA_1, DELTA_2, DELTA_3, it.GetSubmeshIndex());
      }

      {
        const Real DELTA_1 =
          i_r.template interpolateInCell<WAVE_3D_R1>(posCell, pos + d_cachedMaterials[it.GetSubmeshIndex()].d_probeC1) -
          currElastic(WAVE_3D_R1);
        const Real DELTA_2 =
          i_r.template interpolateInCell<WAVE_3D_R2>(posCell, pos + d_cachedMaterials[it.GetSubmeshIndex()].d_probeC2) -
          currElastic(WAVE_3D_R2);
        const Real DELTA_3 =
          i_r.template interpolateInCell<WAVE_3D_R3>(posCell, pos + d_cachedMaterials[it.GetSubmeshIndex()].d_probeC2) -
          currElastic(WAVE_3D_R3);
        nextElastic += incRightWaves(DELTA_1, DELTA_2, DELTA_3, it.GetSubmeshIndex());
      }
    }
  } // for (int i = 0; i < n; i++)
} // BasicGcCalculator3D::compute()


template <typename Space>
template <typename It, typename Reconstructor>
void ContactGcCalculator3D<Space>::computeOuterNodesWaves(It i_begin, It i_end,
  const Reconstructor & i_r) const
{
  const long int n = (long int)(i_end - i_begin);
  #pragma omp parallel for
  for (long int i = 0; i < n; i++)
  {
    It contact = i_begin + i;

    size_t negCell = size_t(-1);
    size_t posCell = size_t(-1);

    size_t contactSides[2];
    contactSides[0] = size_t(-1);
    contactSides[1] = size_t(-1);

    Vector3 contactPoints[2];
    Elastic contactElastics[2];

    for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
    {
      size_t localNegCell, localPosCell;
      contact.findCellAlongDir(d_n, localNegCell, localPosCell, contactSideIndex);

      if(localNegCell != size_t(-1))
      {
        negCell = localNegCell;
        contactSides[0] = contactSideIndex;
        contactPoints[0] = contact.currPos(contactSideIndex);
        contactElastics[0] = contact.currValue(contactSideIndex);
      }
      if(localPosCell != size_t(-1))
      {
        posCell = localPosCell;
        contactSides[1] = contactSideIndex;
        contactPoints[1] = contact.currPos(contactSideIndex);
        contactElastics[1] = contact.currValue(contactSideIndex);
      }
    }
    if(contactSides[0] != size_t(-1) && contact.GetContactTypeIndex(contactSides[0], -1) != size_t(-1))
    {
      if((contact.currNormal(contactSides[0], -1) * (-d_n)) > 0)
      {
        contactSides[0] = size_t(-1);
      }
    }
    if(contactSides[1] != size_t(-1) && contact.GetContactTypeIndex(contactSides[1], -1) != size_t(-1))
    {
      if((contact.currNormal(contactSides[1], -1) *  d_n) > 0)
      {
        contactSides[1] = size_t(-1);
      }
    }

//    if(contact.GetContactSidesCount() == 1 && contactSides[0] != contactSides[1]) printf("wrong cavity boundary");

    if(contactSides[0] == contactSides[1] && contactSides[0] != size_t(-1)) continue; //node counts as internal, should never happen due to boundary hit correction

    for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
    {
      contact.nextValue(contactSideIndex) = contact.tempValue(contactSideIndex);
    }

    if((contactSides[0] != size_t(-1)) && (contactSides[1] != size_t(-1))) //contact
    {
      size_t submeshIndex0 = contact.GetSubmeshIndex(contactSides[0]);
      size_t submeshIndex1 = contact.GetSubmeshIndex(contactSides[1]);

      //normal contact
//      if( (d_cachedMaterials[submeshIndex0].stepMultiplier == d_cachedMaterials[submeshIndex1].stepMultiplier) &&
//          (globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0))
      Real L_DELTA_1 = 0, L_DELTA_2 = 0, L_DELTA_3 = 0;
      Real R_DELTA_1 = 0, R_DELTA_2 = 0, R_DELTA_3 = 0;
      if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0)
      {
        L_DELTA_1 =
          i_r.template interpolateInCell<WAVE_3D_L1>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC1) -
          contactElastics[0](WAVE_3D_L1);
        L_DELTA_2 =
          i_r.template interpolateInCell<WAVE_3D_L2>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2) -
          contactElastics[0](WAVE_3D_L2);
        L_DELTA_3 =
          i_r.template interpolateInCell<WAVE_3D_L3>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2) -
          contactElastics[0](WAVE_3D_L3);
      }else
      if( (d_cachedMaterials[submeshIndex0].stepMultiplier == 
           d_cachedMaterials[submeshIndex1].stepMultiplier * 2) &&
          (globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0))
      {
        L_DELTA_1 =
          i_r.template interpolateInCell<WAVE_3D_L1>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC1 * Real(0.5)) -
          contactElastics[0](WAVE_3D_L1);
        L_DELTA_2 =
          i_r.template interpolateInCell<WAVE_3D_L2>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2 * Real(0.5)) -
          contactElastics[0](WAVE_3D_L2);
        L_DELTA_3 =
          i_r.template interpolateInCell<WAVE_3D_L3>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2 * Real(0.5)) -
          contactElastics[0](WAVE_3D_L3);
      }

      if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0)
      {
        R_DELTA_1 =
          i_r.template interpolateInCell<WAVE_3D_R1>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC1) -
          contactElastics[1](WAVE_3D_R1);
        R_DELTA_2 =
          i_r.template interpolateInCell<WAVE_3D_R2>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2) -
          contactElastics[1](WAVE_3D_R2);
        R_DELTA_3 =
          i_r.template interpolateInCell<WAVE_3D_R3>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2) -
          contactElastics[1](WAVE_3D_R3);
      }else
      if( (d_cachedMaterials[submeshIndex1].stepMultiplier == 
           d_cachedMaterials[submeshIndex0].stepMultiplier * 2) &&
          (globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0))
      {
        R_DELTA_1 =
          i_r.template interpolateInCell<WAVE_3D_R1>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC1 * Real(0.5)) -
          contactElastics[1](WAVE_3D_R1);
        R_DELTA_2 =
          i_r.template interpolateInCell<WAVE_3D_R2>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2 * Real(0.5)) -
          contactElastics[1](WAVE_3D_R2);
        R_DELTA_3 =
          i_r.template interpolateInCell<WAVE_3D_R3>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2 * Real(0.5)) -
          contactElastics[1](WAVE_3D_R3);
      }

      //Elastic groupElastic;
      size_t referenceSide = size_t(-1);
      if(d_cachedMaterials[submeshIndex0].stepMultiplier <
         d_cachedMaterials[submeshIndex1].stepMultiplier)
      {
        if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier != 0)
          continue;
        contact.nextValue(contactSides[0]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
        contact.nextValue(contactSides[0]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);
        if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0)
        {
          contact.nextValue(contactSides[1]) = contact.nextValue(contactSides[0]);
        }
        //referenceSide = contactSides[0];
        //groupElastic = contact.nextValue(contactSides[0]);
      }else
      if(d_cachedMaterials[submeshIndex0].stepMultiplier >
         d_cachedMaterials[submeshIndex1].stepMultiplier)
      {
        if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier != 0)
          continue;
        contact.nextValue(contactSides[1]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
        contact.nextValue(contactSides[1]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);

        if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0)
        {
          contact.nextValue(contactSides[0]) = contact.nextValue(contactSides[1]);
        }
        //referenceSide = contactSides[1];
        //groupElastic = contact.nextValue(contactSides[1]);
      }else
      {
        assert((d_cachedMaterials[submeshIndex0].stepMultiplier ==
                d_cachedMaterials[submeshIndex1].stepMultiplier));
        if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier != 0)
          continue;

        contact.nextValue(contactSides[0]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
        contact.nextValue(contactSides[1]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);

        size_t typeIndex = contact.GetContactTypeIndex(contactSides[0], contactSides[1]);
        if(typeIndex == size_t(-1))
        {
          printf("undefined contact condition\n");
          continue;
        }
        switch(contactSettings[typeIndex].type)
        {
          case ContactSettings::Glue:
            //basic_contact_glue_workspace< space3D<Real> > contactCorrector;
            //contactCorrector.material(d_cachedMaterials[submeshIndex0].d_material, d_cachedMaterials[submeshIndex1].d_material);
            //contactCorrector.correct(contact.nextValue(contactSides[0]), contact.nextValue(contactSides[1]), contact.currNormal(contactSides[0], contactSides[1]));
            contact.nextValue(contactSides[1]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
            contact.nextValue(contactSides[0]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);

            //groupElastic = contact.nextValue(contactSides[0]);
          break;
          case ContactSettings::Slide:
            /*basic_contact_glide_workspace< space3D<Real> > contactCorrector;
            //contactCorrector.material(d_cachedMaterials[submeshIndex0].d_material, d_cachedMaterials[submeshIndex1].d_material);
            basic_border_mixed_workspace< space3D<Real> > space0;
            basic_border_mixed_workspace< space3D<Real> > space1;

            space0.material(d_cachedMaterials[submeshIndex0].d_material);
            space1.material(d_cachedMaterials[submeshIndex1].d_material);
            contactCorrector.correct(contact.nextValue(contactSides[0]), space0, contact.nextValue(contactSides[1]), space1, contact.currNormal(contactSides[0], contactSides[1]));*/
            //contact.nextValue(contactSides[1]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
            //contact.nextValue(contactSides[0]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);

            //groupElastic = contact.nextValue(contactSides[0]);
          break;
        }
      }
      for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
      {
        if((contactSideIndex != contactSides[0]) && (contactSideIndex != contactSides[1]) &&
           (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSideIndex)].stepMultiplier == 0))
        {
          if((contactSettings[contact.GetContactTypeIndex(contactSides[0], contactSideIndex)].type == ContactSettings::Glue) &&
             (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSides[0])].stepMultiplier == 0))
          {
            contact.nextValue(contactSideIndex) = contact.nextValue(contactSides[0]);
          }else
          if((contactSettings[contact.GetContactTypeIndex(contactSides[1], contactSideIndex)].type == ContactSettings::Glue) &&
             (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSides[1])].stepMultiplier == 0))
          {
            contact.nextValue(contactSideIndex) = contact.nextValue(contactSides[1]);
          }
          //contact.nextValue(contactSideIndex) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, contact.GetSubmeshIndex(contactSideIndex));
          //contact.nextValue(contactSideIndex) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, contact.GetSubmeshIndex(contactSideIndex));
          //contact.nextValue(contactSideIndex) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
          //contact.nextValue(contactSideIndex) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);
          //contact.nextValue(contactSideIndex) = groupElastic;
        }
      }
    }else //border or corner
    {
      size_t boundarySide = size_t(-1);

      Real L_DELTA_1 = 0, L_DELTA_2 = 0, L_DELTA_3 = 0;
      Real R_DELTA_1 = 0, R_DELTA_2 = 0, R_DELTA_3 = 0;
      if((contactSides[0] != size_t(-1)) && (contactSides[1] == size_t(-1)))
      {
        size_t submeshIndex0 = contact.GetSubmeshIndex(contactSides[0]);
        if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0)
        {
          L_DELTA_1 =
            i_r.template interpolateInCell<WAVE_3D_L1>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC1) -
            contactElastics[0](WAVE_3D_L1);
          L_DELTA_2 =
            i_r.template interpolateInCell<WAVE_3D_L2>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2) -
            contactElastics[0](WAVE_3D_L2);
          L_DELTA_3 =
            i_r.template interpolateInCell<WAVE_3D_L3>(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2) -
            contactElastics[0](WAVE_3D_L3);
          R_DELTA_1 = 0;//L_DELTA_1;
          R_DELTA_2 = 0;//L_DELTA_2;
          R_DELTA_3 = 0;//L_DELTA_3;

          contact.nextValue(contactSides[0]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
          contact.nextValue(contactSides[0]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex0);

          size_t typeIndex = contact.GetContactTypeIndex(contactSides[0], -1);

          if(typeIndex != size_t(-1))
          {
            switch(boundarySettings[typeIndex].type)
            {
              case BoundarySettings::Free:
              {
                /*basic_border_force_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);*/
                //borderCorrector.correct(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), Vector(0, 0, 0), 0);
                correctFreeBoundary(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), d_cachedMaterials[submeshIndex0].d_material, 0);

                /*ElasticForceCorrector<Space> borderCorrector;
                borderCorrector.setForce(Vector(0, 0, 0));
                borderCorrector.setMaterial(&(d_cachedMaterials[submeshIndex0].d_material));
                borderCorrector.correct(&(contact.nextValue(contactSides[0])), contact.currNormal(contactSides[0], -1));*/
              }
              break;
              case BoundarySettings::Fixed:
              {
                /*basic_border_velocity_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), Vector(0, 0, 0), 0);*/
              }
              break;
              case BoundarySettings::Mixed:
              {
                /*Vector velocity = Vector(0, 0, 0);
                basic_border_mixed_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), velocity, 0, 0);*/
              }
              break;
              case BoundarySettings::Absorb:
              {
                /*basic_border_absorption_workspace< space3D<Real> > borderCorrector;

                basic_border_velocity_workspace< space3D<Real> > velocityCorrector;
                velocityCorrector.material(d_cachedMaterials[submeshIndex0].d_material);

                basic_border_force_workspace< space3D<Real> > forceCorrector;
                forceCorrector.material(d_cachedMaterials[submeshIndex0].d_material);

                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);
                borderCorrector.correct(forceCorrector, velocityCorrector, contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), 0);*/
              }
              break;
              case BoundarySettings::Undefined:
              {
              }
              break;
            }
          }else
          {
            //printf("undefined boundary condition\n");
            continue;
          }
        }

        boundarySide = contactSides[0];
        //groupElastic = contact.nextValue(contactSides[0]);
      }else
      if((contactSides[1] != size_t(-1)) && (contactSides[0] == size_t(-1)))
      {

        size_t submeshIndex1 = contact.GetSubmeshIndex(contactSides[1]);
        if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0)
        {
          size_t submeshIndex1 = contact.GetSubmeshIndex(contactSides[1]);
          R_DELTA_1 =
            i_r.template interpolateInCell<WAVE_3D_R1>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC1) -
            contactElastics[1](WAVE_3D_R1);
          R_DELTA_2 =
            i_r.template interpolateInCell<WAVE_3D_R2>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2) -
            contactElastics[1](WAVE_3D_R2);
          R_DELTA_3 =
            i_r.template interpolateInCell<WAVE_3D_R3>(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2) -
            contactElastics[1](WAVE_3D_R3);
          L_DELTA_1 = 0;//R_DELTA_1;
          L_DELTA_2 = 0;//R_DELTA_2;
          L_DELTA_3 = 0;//R_DELTA_3;

          contact.nextValue(contactSides[1]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex1);
          contact.nextValue(contactSides[1]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);

          size_t typeIndex = contact.GetContactTypeIndex(contactSides[1], -1);

          if(typeIndex != size_t(-1))
          {
            switch(boundarySettings[typeIndex].type)
            {
              case BoundarySettings::Free:
              {
                /*basic_border_force_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);*/
                //borderCorrector.correct(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), Vector(0, 0, 0), 0);

                correctFreeBoundary(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), d_cachedMaterials[submeshIndex1].d_material, 1);

                /*ElasticForceCorrector<Space> borderCorrector;
                borderCorrector.setForce(Vector(0, 0, 0));
                borderCorrector.setMaterial(&(d_cachedMaterials[submeshIndex1].d_material));
                borderCorrector.correct(&(contact.nextValue(contactSides[1])), contact.currNormal(contactSides[1], -1));*/
              }
              break;
              case BoundarySettings::Fixed:
              {
                /*basic_border_velocity_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), Vector(0, 0, 0), 0);*/
              }
              break;
              case BoundarySettings::Mixed:
              {
                /*Vector velocity = Vector(0, 0, 0);
                basic_border_mixed_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), velocity, 0, 0);*/
              }
              case BoundarySettings::Absorb:
              {
                /*basic_border_absorption_workspace< space3D<Real> > borderCorrector;

                basic_border_velocity_workspace< space3D<Real> > velocityCorrector;
                velocityCorrector.material(d_cachedMaterials[submeshIndex1].d_material);

                basic_border_force_workspace< space3D<Real> > forceCorrector;
                forceCorrector.material(d_cachedMaterials[submeshIndex1].d_material);

                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);
                borderCorrector.correct(forceCorrector, velocityCorrector, contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), 0);*/
              }
              break;
              case BoundarySettings::Undefined:
              {
              }
              break;
            }
          }else
          {
            //printf("undefined boundary condition\n");
            continue;
          } 
        }

        boundarySide = contactSides[1];
//        groupElastic = contact.nextValue(contactSides[1]);
      }else
      {
        //corner point
        boundarySide = size_t(-1);
        for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
        {
          size_t submeshIndex = contact.GetSubmeshIndex(contactSideIndex);
          if((globalStepNumber % d_cachedMaterials[submeshIndex].stepMultiplier == 0) &&
             (contact.GetContactTypeIndex(contactSideIndex, -1) != size_t(-1)))
          { 
            boundarySide = contactSideIndex;
            break;
          }
        }

        if(boundarySide == size_t(-1)) continue;

        bool leftBoundary = 0;
        if((d_n * contact.currNormal(boundarySide, -1)) < 0)
        {
          leftBoundary = 1;
        }

        size_t typeIndex = contact.GetContactTypeIndex(boundarySide, -1);

        if(typeIndex != size_t(-1))
        {
          switch(boundarySettings[typeIndex].type)
          {
            case BoundarySettings::Free:
            {
              /*basic_border_force_workspace< space3D<Real> > borderCorrector;
              borderCorrector.material(d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].d_material);*/
              //borderCorrector.correct(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), Vector(0, 0, 0), 0);

              correctFreeBoundary(contact.nextValue(boundarySide), contact.currNormal(boundarySide, -1), d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].d_material, leftBoundary);

              /*ElasticForceCorrector<Space> borderCorrector;
              borderCorrector.setForce(Vector(0, 0, 0));
              borderCorrector.setMaterial(&(d_cachedMaterials[submeshIndex1].d_material));
              borderCorrector.correct(&(contact.nextValue(contactSides[1])), contact.currNormal(contactSides[1], -1));*/
            }
            break;
            case BoundarySettings::Fixed:
            {
              /*basic_border_velocity_workspace< space3D<Real> > borderCorrector;
              borderCorrector.material(d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].d_material);
              borderCorrector.correct(contact.nextValue(boundarySide), contact.currNormal(boundarySide, -1), Vector(0, 0, 0), 0);*/
            }
            break;
            case BoundarySettings::Mixed:
            {
              /*Vector velocity = Vector(0, 0, 0);
              basic_border_mixed_workspace< space3D<Real> > borderCorrector;
              borderCorrector.material(d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].d_material);
              borderCorrector.correct(contact.nextValue(boundarySide), contact.currNormal(boundarySide, -1), velocity, 0, 0);*/
            }
            case BoundarySettings::Absorb:
            {
              /*basic_border_absorption_workspace< space3D<Real> > borderCorrector;

              basic_border_velocity_workspace< space3D<Real> > velocityCorrector;
              velocityCorrector.material(d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].d_material);

              basic_border_force_workspace< space3D<Real> > forceCorrector;
              forceCorrector.material(d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].d_material);

              borderCorrector.material(d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].d_material);
              borderCorrector.correct(forceCorrector, velocityCorrector, contact.nextValue(boundarySide), contact.currNormal(boundarySide, -1), 0);*/
            }
            break;
            case BoundarySettings::Undefined:
            {
            }
            break;
          }
        }
      }
      if(boundarySide == size_t(-1)) continue; //no appropriate boundary to copy from

      for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
      {
        if(globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSideIndex)].stepMultiplier == 0)
        {
          if((boundarySide != contactSideIndex) &&
             (contactSettings[contact.GetContactTypeIndex(boundarySide, contactSideIndex)].type == ContactSettings::Glue) &&
             (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].stepMultiplier == 0))
          {
            contact.nextValue(contactSideIndex) = contact.nextValue(boundarySide);
          }
          //contact.nextValue(contactSideIndex) = groupElastic;
        }
      }
       //border
    }
  } // for (int i = 0; i < n; i++)
} // BasicGcCalculator3D::computeOuterNodes()

template <typename Space>
typename ContactGcCalculator3D<Space>::Elastic
  ContactGcCalculator3D<Space>::incDirect(const Elastic &lc1, const Elastic &lc2, const Elastic &rc1, const Elastic &rc2, const size_t i_materialId) const
{
  return incLeftDirect(lc1, lc2, i_materialId) + incRightDirect(rc1, rc2, i_materialId);
}

template <typename Space>
typename ContactGcCalculator3D<Space>::Elastic
  ContactGcCalculator3D<Space>::incLeftDirect(const Elastic &lc1, const Elastic &lc2, const size_t i_materialId) const
{
  Elastic res;

  res.SetVelocity(Vector3(0, 0, 0));
  res.SetTension(Tensor3(0));

  const CachedMaterial &material = d_cachedMaterials[i_materialId];

  //+-c1
  {
    Real leftMult1  = Real(0.5) * ((d_n * lc1.GetVelocity()) - DoubleConvolution(d_N, lc1.GetTension()) * material.d_rrhoc1);
    res.SetVelocity(res.GetVelocity() + leftMult1  * d_n);
    res.SetTension (res.GetTension () + Tensor3(leftMult1 * material.d_rhoc3));
    res.SetTension (res.GetTension () + d_N * (-leftMult1  * material.d_rhoc1mc3));
  }

  //+-c2
  {
    Real leftMult1 = Real(0.5) * (d_n1 * lc2.GetVelocity() - DoubleConvolution(d_N1, lc2.GetTension()) * material.d_rrhoc2);
    res.SetVelocity(res.GetVelocity() + leftMult1 * d_n1);
    res.SetTension (res.GetTension () + d_N1 * Real(-2.0) * leftMult1 * material.d_rhoc2);

    Real leftMult2 = Real(0.5) * (d_n2 * lc2.GetVelocity() - DoubleConvolution(d_N2, lc2.GetTension()) * material.d_rrhoc2);
    res.SetVelocity(res.GetVelocity() + leftMult2 * d_n2);
    res.SetTension (res.GetTension () + d_N2 * Real(-2.0) * leftMult2 * material.d_rhoc2);
  }
  return res;
}

template <typename Space>
typename ContactGcCalculator3D<Space>::Elastic
  ContactGcCalculator3D<Space>::incRightDirect(const Elastic &rc1, const Elastic &rc2, const size_t i_materialId) const
{
  Elastic res;

  res.SetVelocity(Vector3(0, 0, 0));
  res.SetTension(Tensor3(0));

  const CachedMaterial &material = d_cachedMaterials[i_materialId];

  //+-c1
  {
    Real rightMult1 = Real(0.5) * (d_n * rc1.GetVelocity() + DoubleConvolution(d_N, rc1.GetTension()) * material.d_rrhoc1);
    res.SetVelocity(res.GetVelocity() + rightMult1 * d_n);
    res.SetTension (res.GetTension () + Tensor3(rightMult1 * material.d_rhoc3));
    res.SetTension (res.GetTension () + d_N * rightMult1 * material.d_rhoc1mc3);
  }

  //+-c2
  {
    Real rightMult1 = Real(0.5) * (d_n1 * rc2.GetVelocity() + DoubleConvolution(d_N1, rc2.GetTension()) * material.d_rrhoc2);
    res.SetVelocity(res.GetVelocity() + rightMult1 * d_n1);
    res.SetTension (res.GetTension () + d_N1 * Real(2.0) * rightMult1 * material.d_rhoc2);

    Real rightMult2 = Real(0.5) * (d_n2 * rc2.GetVelocity() + DoubleConvolution(d_N2, rc2.GetTension()) * material.d_rrhoc2);
    res.SetVelocity(res.GetVelocity() + rightMult2 * d_n2);
    res.SetTension (res.GetTension () + d_N2 * Real(2.0) * rightMult2 * material.d_rhoc2);
  }
  return res;
}


template <typename Space>
template <typename It, typename Reconstructor>
void ContactGcCalculator3D<Space>::computeInnerNodesDirect(It i_begin, It i_end,
  const Reconstructor & i_r) const
{
  const long int n = (long int)(i_end - i_begin);
  #pragma omp parallel for
  for (long int i = 0; i < n; i++)
  {
    It it = i_begin + i;
    const Elastic & currElastic = it.currValue();
    Elastic & nextElastic = it.nextValue();
    Vector3 pos = it.currPos();



    size_t negCell, posCell;
    it.findCellAlongDir(d_n, negCell, posCell);

    if ((negCell != (size_t)(-1)) && (posCell != (size_t)(-1)))
    {
      Elastic lc1 = i_r.interpolateElastic(negCell, pos - d_cachedMaterials[it.GetSubmeshIndex()].d_probeC1) - it.currValue();
      Elastic lc2 = i_r.interpolateElastic(negCell, pos - d_cachedMaterials[it.GetSubmeshIndex()].d_probeC2) - it.currValue();
      Elastic rc1 = i_r.interpolateElastic(posCell, pos + d_cachedMaterials[it.GetSubmeshIndex()].d_probeC1) - it.currValue();
      Elastic rc2 = i_r.interpolateElastic(posCell, pos + d_cachedMaterials[it.GetSubmeshIndex()].d_probeC2) - it.currValue();
      nextElastic += incDirect(lc1, lc2, rc1, rc2, it.GetSubmeshIndex());
    }
  } // for (int i = 0; i < n; i++)
}


template <typename Space>
template <typename It, typename Reconstructor>
void ContactGcCalculator3D<Space>::computeOuterNodesDirect(It i_begin, It i_end,
  const Reconstructor & i_r) const
{
  const long int n = (long int)(i_end - i_begin);
  #pragma omp parallel for
  for (long int i = 0; i < n; i++)
  {
    It contact = i_begin + i;

    size_t negCell = size_t(-1);
    size_t posCell = size_t(-1);

    size_t contactSides[2];
    contactSides[0] = size_t(-1);
    contactSides[1] = size_t(-1);

    Vector3 contactPoints[2];
    Elastic contactElastics[2];

    for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
    {
      size_t localNegCell, localPosCell;
      contact.findCellAlongDir(d_n, localNegCell, localPosCell, contactSideIndex);

      if(localNegCell != size_t(-1))
      {
        negCell = localNegCell;
        contactSides[0] = contactSideIndex;
        contactPoints[0] = contact.currPos(contactSideIndex);
        contactElastics[0] = contact.currValue(contactSideIndex);
      }
      if(localPosCell != size_t(-1))
      {
        posCell = localPosCell;
        contactSides[1] = contactSideIndex;
        contactPoints[1] = contact.currPos(contactSideIndex);
        contactElastics[1] = contact.currValue(contactSideIndex);
      }
    }
    /*if(contactSides[0] != size_t(-1) && contact.GetContactTypeIndex(contactSides[0], -1) != size_t(-1))
    {
      if(scalar_product(contact.currNormal(contactSides[0], -1), -d_n) > 0)
      {
        contactSides[0] = size_t(-1);
      }
    }
    if(contactSides[1] != size_t(-1) && contact.GetContactTypeIndex(contactSides[1], -1) != size_t(-1))
    {
      if(scalar_product(contact.currNormal(contactSides[1], -1),  d_n) > 0)
      {
        contactSides[1] = size_t(-1);
      }
    }*/

    if(contactSides[0] == contactSides[1]) continue;

    if((contactSides[0] != size_t(-1)) && (contactSides[1] != size_t(-1))) //contact
    {
      size_t submeshIndex0 = contact.GetSubmeshIndex(contactSides[0]);
      size_t submeshIndex1 = contact.GetSubmeshIndex(contactSides[1]);

      if(contactSides[0] != contactSides[1]) //point is not internal
      {
        //normal contact
  //      if( (d_cachedMaterials[submeshIndex0].stepMultiplier == d_cachedMaterials[submeshIndex1].stepMultiplier) &&
  //          (globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0))
        Elastic lc1, lc2, rc1, rc2;
        lc1.SetVelocity(Vector3(0, 0, 0));
        lc2.SetVelocity(Vector3(0, 0, 0));
        rc1.SetVelocity(Vector3(0, 0, 0));
        rc2.SetVelocity(Vector3(0, 0, 0));

        lc1.SetTension(Tensor3(0));
        lc2.SetTension(Tensor3(0));
        rc1.SetTension(Tensor3(0));
        rc2.SetTension(Tensor3(0));

        if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0)
        {
          lc1 = i_r.interpolateElastic(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC1) - contactElastics[0];
          lc2 = i_r.interpolateElastic(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2) - contactElastics[0];
        }else
        if( (d_cachedMaterials[submeshIndex0].stepMultiplier == 
             d_cachedMaterials[submeshIndex1].stepMultiplier * 2) &&
            (globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0))
        {
          lc1 = i_r.interpolateElastic(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC1 * Real(0.5)) - contactElastics[0];
          lc2 = i_r.interpolateElastic(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2 * Real(0.5)) - contactElastics[0];
        }

        if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0)
        {
          rc1 = i_r.interpolateElastic(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC1) - contactElastics[1];
          rc2 = i_r.interpolateElastic(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2) - contactElastics[1];
        }else
        if( (d_cachedMaterials[submeshIndex1].stepMultiplier == 
             d_cachedMaterials[submeshIndex0].stepMultiplier * 2) &&
            (globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0))
        {
          rc1 = i_r.interpolateElastic(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC1 * Real(0.5)) - contactElastics[1];
          rc2 = i_r.interpolateElastic(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2 * Real(0.5)) - contactElastics[1];
        }

        //Elastic groupElastic;
        size_t referenceSide = size_t(-1);
        if(d_cachedMaterials[submeshIndex0].stepMultiplier <
           d_cachedMaterials[submeshIndex1].stepMultiplier)
        {
          if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier != 0)
            continue;
          contact.nextValue(contactSides[0]) += incLeftDirect (lc1, lc2, submeshIndex0);
          contact.nextValue(contactSides[0]) += incRightDirect(rc1, rc2, submeshIndex1);
          if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0)
          {
            contact.nextValue(contactSides[1]) = contact.nextValue(contactSides[0]);
          }
          //referenceSide = contactSides[0];
          //groupElastic = contact.nextValue(contactSides[0]);
        }else
        if(d_cachedMaterials[submeshIndex0].stepMultiplier >
           d_cachedMaterials[submeshIndex1].stepMultiplier)
        {
          if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier != 0)
            continue;
          contact.nextValue(contactSides[1]) += incLeftDirect (lc1, lc2, submeshIndex0);
          contact.nextValue(contactSides[1]) += incRightDirect(rc1, rc2, submeshIndex1);

          if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0)
          {
            contact.nextValue(contactSides[0]) = contact.nextValue(contactSides[1]);
          }
          //referenceSide = contactSides[1];
          //groupElastic = contact.nextValue(contactSides[1]);
        }else
        {
          assert((d_cachedMaterials[submeshIndex0].stepMultiplier ==
                  d_cachedMaterials[submeshIndex1].stepMultiplier));
          if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier != 0)
            continue;

          contact.nextValue(contactSides[0]) += incLeftDirect (lc1, lc2, submeshIndex0);
          contact.nextValue(contactSides[1]) += incRightDirect(rc1, rc2, submeshIndex1);

          size_t typeIndex = contact.GetContactTypeIndex(contactSides[0], contactSides[1]);
          if(typeIndex == size_t(-1))
          {
            printf("undefined contact condition\n");
            continue;
          }
          switch(contactSettings[typeIndex].type)
          {
            case ContactSettings::Glue:
              //basic_contact_glue_workspace< space3D<Real> > contactCorrector;
              //contactCorrector.material(d_cachedMaterials[submeshIndex0].d_material, d_cachedMaterials[submeshIndex1].d_material);
              //contactCorrector.correct(contact.nextValue(contactSides[0]), contact.nextValue(contactSides[1]), contact.currNormal(contactSides[0], contactSides[1]));
              contact.nextValue(contactSides[1]) += incLeftDirect (lc1, lc2, submeshIndex0);
              contact.nextValue(contactSides[0]) += incRightDirect(rc1, rc2, submeshIndex1);

              //groupElastic = contact.nextValue(contactSides[0]);
            break;
            case ContactSettings::Slide:
              /*basic_contact_glide_workspace< space3D<Real> > contactCorrector;
              //contactCorrector.material(d_cachedMaterials[submeshIndex0].d_material, d_cachedMaterials[submeshIndex1].d_material);
              basic_border_mixed_workspace< space3D<Real> > space0;
              basic_border_mixed_workspace< space3D<Real> > space1;

              space0.material(d_cachedMaterials[submeshIndex0].d_material);
              space1.material(d_cachedMaterials[submeshIndex1].d_material);
              contactCorrector.correct(contact.nextValue(contactSides[0]), space0, contact.nextValue(contactSides[1]), space1, contact.currNormal(contactSides[0], contactSides[1]));*/
              //contact.nextValue(contactSides[1]) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
              //contact.nextValue(contactSides[0]) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);

              //groupElastic = contact.nextValue(contactSides[0]);
            break;
          }
        }
      }
      for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
      {
        if((contactSideIndex != contactSides[0]) && (contactSideIndex != contactSides[1]) &&
           (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSideIndex)].stepMultiplier == 0))
        {
          if((contactSettings[contact.GetContactTypeIndex(contactSides[0], contactSideIndex)].type == ContactSettings::Glue) &&
             (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSides[0])].stepMultiplier == 0))
          {
            contact.nextValue(contactSideIndex) = contact.nextValue(contactSides[0]);
          }else
          if((contactSettings[contact.GetContactTypeIndex(contactSides[1], contactSideIndex)].type == ContactSettings::Glue) &&
             (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSides[1])].stepMultiplier == 0))
          {
            contact.nextValue(contactSideIndex) = contact.nextValue(contactSides[1]);
          }
          //contact.nextValue(contactSideIndex) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, contact.GetSubmeshIndex(contactSideIndex));
          //contact.nextValue(contactSideIndex) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, contact.GetSubmeshIndex(contactSideIndex));
          //contact.nextValue(contactSideIndex) += incLeftWaves (L_DELTA_1, L_DELTA_2, L_DELTA_3, submeshIndex0);
          //contact.nextValue(contactSideIndex) += incRightWaves(R_DELTA_1, R_DELTA_2, R_DELTA_3, submeshIndex1);
          //contact.nextValue(contactSideIndex) = groupElastic;
        }
      }
    }else //border or corner
    {
      Elastic groupElastic;
      size_t boundarySide = size_t(-1);

      Elastic lc1, lc2, rc1, rc2;
      lc1.SetVelocity(Vector3(0, 0, 0));
      lc2.SetVelocity(Vector3(0, 0, 0));
      rc1.SetVelocity(Vector3(0, 0, 0));
      rc2.SetVelocity(Vector3(0, 0, 0));

      lc1.SetTension(Tensor3(0));
      lc2.SetTension(Tensor3(0));
      rc1.SetTension(Tensor3(0));
      rc2.SetTension(Tensor3(0));
      if((contactSides[0] != size_t(-1)) && (contactSides[1] == size_t(-1)))
      {
        size_t submeshIndex0 = contact.GetSubmeshIndex(contactSides[0]);
        if(globalStepNumber % d_cachedMaterials[submeshIndex0].stepMultiplier == 0)
        {
          lc1 =
            i_r.interpolateElastic(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC1) -
            contactElastics[0];
          lc2 =
            i_r.interpolateElastic(negCell, contactPoints[0] - d_cachedMaterials[submeshIndex0].d_probeC2) -
            contactElastics[0];

          contact.nextValue(contactSides[0]) += incLeftDirect (lc1, lc2, submeshIndex0);

          size_t typeIndex = contact.GetContactTypeIndex(contactSides[0], -1);

          if(typeIndex != size_t(-1))
          {
            switch(boundarySettings[typeIndex].type)
            {
              case BoundarySettings::Free:
              {
                /*basic_border_force_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), Vector(0, 0, 0), 0);*/
                correctFreeBoundary(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), d_cachedMaterials[submeshIndex0].d_material, 0);

                /*ElasticForceCorrector<Space> borderCorrector;
                borderCorrector.setForce(Vector(0, 0, 0));
                borderCorrector.setMaterial(&(d_cachedMaterials[submeshIndex0].d_material));
                borderCorrector.correct(&(contact.nextValue(contactSides[0])), contact.currNormal(contactSides[0], -1));*/
              }
              break;
              case BoundarySettings::Fixed:
              {
                /*basic_border_velocity_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), Vector(0, 0, 0), 0);*/
              }
              break;
              case BoundarySettings::Mixed:
              {
                /*Vector velocity = Vector(0, 0, 0);
                basic_border_mixed_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), velocity, 0, 0);*/
              }
              break;
              case BoundarySettings::Absorb:
              {
                /*basic_border_absorption_workspace< space3D<Real> > borderCorrector;

                basic_border_velocity_workspace< space3D<Real> > velocityCorrector;
                velocityCorrector.material(d_cachedMaterials[submeshIndex0].d_material);

                basic_border_force_workspace< space3D<Real> > forceCorrector;
                forceCorrector.material(d_cachedMaterials[submeshIndex0].d_material);

                borderCorrector.material(d_cachedMaterials[submeshIndex0].d_material);
                borderCorrector.correct(forceCorrector, velocityCorrector, contact.nextValue(contactSides[0]), contact.currNormal(contactSides[0], -1), 0);*/
              }
              break;
              case BoundarySettings::Undefined:
              {
              }
              break;
            }
          }else
          {
            //printf("undefined boundary condition\n");
            continue;
          }
        }

        boundarySide = contactSides[0];
        //groupElastic = contact.nextValue(contactSides[0]);
      }else
      if((contactSides[1] != size_t(-1)) && (contactSides[0] == size_t(-1)))
      {

        size_t submeshIndex1 = contact.GetSubmeshIndex(contactSides[1]);
        if(globalStepNumber % d_cachedMaterials[submeshIndex1].stepMultiplier == 0)
        {
          size_t submeshIndex1 = contact.GetSubmeshIndex(contactSides[1]);
          rc1 =
            i_r.interpolateElastic(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC1) - contactElastics[1];
          rc2 =
            i_r.interpolateElastic(posCell, contactPoints[1] + d_cachedMaterials[submeshIndex1].d_probeC2) - contactElastics[1];

          contact.nextValue(contactSides[1]) += incRightDirect(rc1, rc2, submeshIndex1);

          size_t typeIndex = contact.GetContactTypeIndex(contactSides[1], -1);

          if(typeIndex != size_t(-1))
          {
            switch(boundarySettings[typeIndex].type)
            {
              case BoundarySettings::Free:
              {
                /*basic_border_force_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), Vector(0, 0, 0), 0);*/

                correctFreeBoundary(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), d_cachedMaterials[submeshIndex1].d_material, 1);

                /*ElasticForceCorrector<Space> borderCorrector;
                borderCorrector.setForce(Vector(0, 0, 0));
                borderCorrector.setMaterial(&(d_cachedMaterials[submeshIndex1].d_material));
                borderCorrector.correct(&(contact.nextValue(contactSides[1])), contact.currNormal(contactSides[1], -1));*/
              }
              break;
              case BoundarySettings::Fixed:
              {
                /*basic_border_velocity_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), Vector(0, 0, 0), 0);*/
              }
              break;
              case BoundarySettings::Mixed:
              {
                /*Vector velocity = Vector(0, 0, 0);
                basic_border_mixed_workspace< space3D<Real> > borderCorrector;
                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);
                borderCorrector.correct(contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), velocity, 0, 0);*/
              }
              case BoundarySettings::Absorb:
              {
                /*basic_border_absorption_workspace< space3D<Real> > borderCorrector;

                basic_border_velocity_workspace< space3D<Real> > velocityCorrector;
                velocityCorrector.material(d_cachedMaterials[submeshIndex1].d_material);

                basic_border_force_workspace< space3D<Real> > forceCorrector;
                forceCorrector.material(d_cachedMaterials[submeshIndex1].d_material);

                borderCorrector.material(d_cachedMaterials[submeshIndex1].d_material);
                borderCorrector.correct(forceCorrector, velocityCorrector, contact.nextValue(contactSides[1]), contact.currNormal(contactSides[1], -1), 0);*/
              }
              break;
              case BoundarySettings::Undefined:
              {
              }
              break;
            }
          }else
          {
            //printf("undefined boundary condition\n");
            continue;
          } 
        }

        boundarySide = contactSides[1];
//        groupElastic = contact.nextValue(contactSides[1]);
      }else
      {
        //corner point

        boundarySide = contactSides[0];
        groupElastic = contact.nextValue(0);
        continue;
      }

      for(size_t contactSideIndex = 0; contactSideIndex < contact.GetContactSidesCount(); contactSideIndex++)
      {
        if(globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(contactSideIndex)].stepMultiplier == 0)
        {
          if((boundarySide != contactSideIndex) &&
             (contactSettings[contact.GetContactTypeIndex(boundarySide, contactSideIndex)].type == ContactSettings::Glue) &&
             (globalStepNumber % d_cachedMaterials[contact.GetSubmeshIndex(boundarySide)].stepMultiplier == 0))
          {
            contact.nextValue(contactSideIndex) = contact.nextValue(boundarySide);
          }
          //contact.nextValue(contactSideIndex) = groupElastic;
        }
      }
    }
  } // for (int i = 0; i < n; i++)
} // BasicGcCalculator3D::computeOuterNodes()


template <typename Space>
template <typename It>
void ContactGcCalculator3D<Space>::copyNextToCurr(It i_begin, It i_end) const
{
  const long int n = (long int)(i_end - i_begin);
  #pragma omp parallel for
  for (long int i = 0; i < n; ++i)
  {
    It it = i_begin + i;
    Elastic & currElastic = it.currValue();
    const Elastic & nextElastic = it.nextValue();
    currElastic = nextElastic;

  }
}



template <typename SpaceT>
typename ContactGcCalculator3D<SpaceT>::Vector3 ContactGcCalculator3D<SpaceT>::getSplitDirection()
{
  return d_n;
}

template <typename Space>
void ContactGcCalculator3D<Space>::generateSplitVectors()
{
  if(!d_randomizeSplitAxes)
  {
    basisVectors[0] = Vector3::xAxis();
    basisVectors[1] = Vector3::yAxis();
    basisVectors[2] = Vector3::zAxis();
  }else
  {
    basisVectors[0] = Vector3(Real(rand()) / Real(RAND_MAX) - Real(0.5), Real(rand()) / Real(RAND_MAX) - Real(0.5), Real(rand()) / Real(RAND_MAX) - Real(0.5));
    if((basisVectors[0]).Len() > Real(0.0))
      basisVectors[0] = basisVectors[0].GetNorm();
    else
      basisVectors[0] = Vector3::xAxis();

    basisVectors[1] = Vector3(Real(rand()) / Real(RAND_MAX) - Real(0.5), Real(rand()) / Real(RAND_MAX) - Real(0.5), Real(rand()) / Real(RAND_MAX) - Real(0.5));
    if((basisVectors[1]).Len() > Real(0.0))
      basisVectors[1] = basisVectors[1].GetNorm();
    else
      basisVectors[1] = Vector3::yAxis();

    basisVectors[2] = (basisVectors[0] ^ basisVectors[1]);
    if((basisVectors[2]).Len() > Real(0.0))
      basisVectors[2] = basisVectors[2].GetNorm();
    else
    {
      basisVectors[0] = Vector3::xAxis();
      basisVectors[1] = Vector3::yAxis();
      basisVectors[2] = Vector3::zAxis();
    }

    basisVectors[1] = basisVectors[0] ^ basisVectors[2];
  }
}


template <typename Space>
void ContactGcCalculator3D<Space>::initiateSplitStep(size_t i_splitPhase, Real i_step, size_t i_globalStepNumber)
{
  assert(i_splitPhase == 0 || i_splitPhase == 1 || i_splitPhase == 2);

  globalStepNumber = i_globalStepNumber;

  d_n   = basisVectors[(i_splitPhase + 0) % 3];
  d_n1  = basisVectors[(i_splitPhase + 1) % 3];
  d_n2  = basisVectors[(i_splitPhase + 2) % 3];

  // calculate cached values to speed up calculations
  d_N = TensorProduct(d_n);
  d_N1 = SymmetricTensorProduct(d_n, d_n1) * Real(0.5);
  d_N2 = SymmetricTensorProduct(d_n, d_n2) * Real(0.5);

  for(size_t materialIndex = 0; materialIndex < d_cachedMaterials.size(); materialIndex++)
  {
    d_cachedMaterials[materialIndex].d_probeC1 =
      d_cachedMaterials[materialIndex].d_material.c1() * i_step * d_n * Real(d_cachedMaterials[materialIndex].stepMultiplier);
    d_cachedMaterials[materialIndex].d_probeC2 =
      d_cachedMaterials[materialIndex].d_material.c2() * i_step * d_n * Real(d_cachedMaterials[materialIndex].stepMultiplier);
  }
}

template <typename Space>
void ContactGcCalculator3D<Space>::SetBoundarySettings(const BoundarySettings boundarySettings, const IndexType typeIndex)
{
  if(this->boundarySettings.size() < typeIndex + 1) this->boundarySettings.resize(typeIndex + 1);
  this->boundarySettings[typeIndex] = boundarySettings;
}

template <typename Space>
void ContactGcCalculator3D<Space>::SetContactSettings (const ContactSettings contactSettings, const IndexType typeIndex)
{
  if(this->contactSettings.size() < typeIndex + 1) this->contactSettings.resize(typeIndex + 1);
  this->contactSettings[typeIndex] = contactSettings;
}


template <typename Space>
void ContactGcCalculator3D<Space>::computeWaves(const Elastic & i_e, const size_t i_materialId,
  Real o_waves[WAVE_3D_COUNT]) const
{
  const Real V_N   = i_e.GetVelocity() * d_n;
  const Real S_N   = d_cachedMaterials[i_materialId].d_rrhoc1 * DoubleConvolution(i_e.GetTension(), d_N);

  const Real V_N_1 = i_e.GetVelocity() * d_n1;
  const Real S_N_1 = d_cachedMaterials[i_materialId].d_rrhoc2 * DoubleConvolution(i_e.GetTension(), d_N1);

  const Real V_N_2 = i_e.GetVelocity() * d_n2;
  const Real S_N_2 = d_cachedMaterials[i_materialId].d_rrhoc2 * DoubleConvolution(i_e.GetTension(), d_N2);

  o_waves[WAVE_3D_L1] = V_N - S_N;
  o_waves[WAVE_3D_R1] = V_N + S_N;
  o_waves[WAVE_3D_L2] = V_N_1 - S_N_1;
  o_waves[WAVE_3D_R2] = V_N_1 + S_N_1;
  o_waves[WAVE_3D_L3] = V_N_2 - S_N_2;
  o_waves[WAVE_3D_R3] = V_N_2 + S_N_2;
}

template <typename Space>
typename ContactGcCalculator3D<Space>::Elastic
ContactGcCalculator3D<Space>::incLeftWaves(
  Real i_delta1, Real i_delta2, Real i_delta3, const size_t i_materialId) const
{
  Tensor3 rhoc3I(d_cachedMaterials[i_materialId].d_rhoc3);
  return Elastic(
    (i_delta1 * d_n + i_delta2 * d_n1 + i_delta3 * d_n2) * Real(0.5),
    (d_N * d_cachedMaterials[i_materialId].d_rhoc1mc3 + rhoc3I) * Real(-0.5) * i_delta1 -
     d_N1 * i_delta2 * d_cachedMaterials[i_materialId].d_rhoc2 - d_N2 * i_delta3 * d_cachedMaterials[i_materialId].d_rhoc2
  );
}

template <typename Space>
typename ContactGcCalculator3D<Space>::Elastic
ContactGcCalculator3D<Space>::incRightWaves(
  Real i_delta1, Real i_delta2, Real i_delta3, const size_t i_materialId) const
{
  Tensor3 rhoc3I(d_cachedMaterials[i_materialId].d_rhoc3);
  return Elastic(
    Real(0.5) * (i_delta1 * d_n + i_delta2 * d_n1 + i_delta3 * d_n2),
    (d_N * d_cachedMaterials[i_materialId].d_rhoc1mc3 + rhoc3I) * Real(0.5) * i_delta1 +
     d_N1 * i_delta2 * d_cachedMaterials[i_materialId].d_rhoc2 + d_N2 * i_delta3 * d_cachedMaterials[i_materialId].d_rhoc2
  );
}

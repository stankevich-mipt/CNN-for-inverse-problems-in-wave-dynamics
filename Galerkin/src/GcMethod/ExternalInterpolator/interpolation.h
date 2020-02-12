#ifndef InterpolationAleanera2012
#define InterpolationAleanera2012
#include "interpolationHeader.h"

///////////////////////////Монотонная/////////////////////

template <typename ValueType>
ValueType interpolate_Compare_Middle(
const int a, const int b, const int c,
const ValueType vA, const ValueType vB,
const char Axis,
const ValueType * R,
const ValueType * V)
//R из 15 точек x, y, z (просто так), xA, yA, zA, xB, yB, zB, xC, ...
//V как в интерполяторе для N=2
//abc(4-a-b-c) идентифицирует номер точки на середние ребра
//vA, vB - то, с чем ее будем сравнивать
{
  ValueType vM;
  ValueType R0 [15];
  int i;

  for(i=0; (i<3); i++)
  {
    R0[i]=R[12+i]+(R[3+i]-R[12+i])*a/(ValueType)4.0+(R[6+i]-R[12+i])*b/(ValueType)4.0+(R[9+i]-R[12+i])*c/(ValueType)4.0;
  };
  for(i=3; (i<15); i++) //i=3
  {
    R0[i]=R[i];
  };

  vM=Interpolate_3D<ValueType>(2, 'P', Axis, R0, V);

  if (vA>vB)
  {
    if ((vM>vA)||(vB>vM))
    {
      return((vA+vB)/2);
    };
  }
  else 
  {
    if ((vM<vA)||(vB<vM))
    {
      return((vA+vB)/2);
    };
  };
  return(vM);
}

template <typename ValueType>
ValueType Interpolate_3D_Mono(
const char Axis,
const ValueType * R,
const ValueType * V)
{
  ValueType V0[35];

/*v4000=V0[0]; 
  v0400=V0[1]; 
  v0040=V0[2]; 
  v0004=V0[3];

  v3100=V0[4];
  v2200=V0[5];
  v1300=V0[6];

  v3010=V0[7];
  v2020=V0[8];
  v1030=V0[9];

  v3001=V0[10];
  v2002=V0[11];
  v1003=V0[12];

  v0310=V0[13]; 
  v0220=V0[14];
  v0130=V0[15];

  v0301=V0[16];
  v0202=V0[17];
  v0103=V0[18];

  v0031=V0[19]; 
  v0022=V0[20]; 
  v0013=V0[21]; 

  v2110=V0[22];
  v1210=V0[23];
  v1120=V0[24];

  v2101=V0[25];
  v1201=V0[26];
  v1102=V0[27];

  v2011=V0[28];
  v1021=V0[29];
  v1012=V0[30];

  v0211=V0[31]; 
  v0121=V0[32]; 
  v0112=V0[33]; 

  v1111=V0[34]; 

  v2000=V[0]; 
  v0200=V[1]; 
  v0020=V[2]; 
  v0002=V[3];

  v1100=V[4];
  v1010=V[5];
  v1001=V[6];
  v0110=V[7];
  v0101=V[8];
  v0011=V[9]; */

  
  V0[ 0]=V[0]; // 4000 2000
  V0[ 1]=V[1]; // 0400 0200
  V0[ 2]=V[2]; // 0040 0020 
  V0[ 3]=V[3]; // 0004 0002

  V0[ 5]=V[4]; // 2200 1100
  V0[ 8]=V[5]; // 2020 1010
  V0[11]=V[6]; // 2002 1001
  V0[14]=V[7]; // 0220 0110
  V0[17]=V[8]; // 0202 0101
  V0[20]=V[9]; // 0022 0011

//  3100 2000 1100  v3100=V0[ 4];  v2000=V[0];  v1100=V[4];
  V0[ 4]=interpolate_Compare_Middle<ValueType>(3, 1, 0, V[0], V[4], Axis, R, V);

//  0310 0200 0110  v0310=V0[13];  v0200=V[1];  v0110=V[7];
  V0[13]=interpolate_Compare_Middle<ValueType>(0, 3, 1, V[1], V[7], Axis, R, V);

//  0031 0020 0011  v0031=V0[19];  v0020=V[2];  v0011=V[9];
  V0[19]=interpolate_Compare_Middle<ValueType>(0, 0, 3, V[2], V[9], Axis, R, V);

//  1003 0002 1001  v1003=V0[12];  v0002=V[3];  v1001=V[6];
  V0[12]=interpolate_Compare_Middle<ValueType>(1, 0, 0, V[3], V[6], Axis, R, V);


//  1300 0200 1100  v1300=V0[ 6];  v0200=V[1];  v1100=V[4];
  V0[ 6]=interpolate_Compare_Middle<ValueType>(1, 3, 0, V[1], V[4], Axis, R, V);

//  0130 0020 0110  v0130=V0[15];  v0020=V[2];  v0110=V[7];
  V0[15]=interpolate_Compare_Middle<ValueType>(0, 1, 3, V[2], V[7], Axis, R, V);

//  0013 0002 0011  v0013=V0[21];  v0002=V[3];  v0011=V[9];
  V0[21]=interpolate_Compare_Middle<ValueType>(0, 0, 1, V[3], V[9], Axis, R, V);

//  3001 2000 1001  v3001=V0[10];  v2000=V[0];  v1001=V[6];
  V0[10]=interpolate_Compare_Middle<ValueType>(3, 0, 0, V[0], V[6], Axis, R, V);


//  3010 2000 1010  v3010=V0[ 7];  v2000=V[0];  v1010=V[5];
  V0[ 7]=interpolate_Compare_Middle<ValueType>(3, 0, 1, V[0], V[5], Axis, R, V);

//  0301 0200 0101  v0301=V0[16];  v0200=V[1];  v0101=V[8];
  V0[16]=interpolate_Compare_Middle<ValueType>(0, 3, 0, V[1], V[8], Axis, R, V);

//  1030 0020 1010  v1030=V0[ 9];  v0020=V[2];  v1010=V[5];
  V0[ 9]=interpolate_Compare_Middle<ValueType>(1, 0, 3, V[2], V[5], Axis, R, V);

//  0103 0002 0101  v0103=V0[18];  v0002=V[3];  v0101=V[8];
  V0[18]=interpolate_Compare_Middle<ValueType>(0, 1, 0, V[3], V[8], Axis, R, V);


//  2110 1100 1010  v2110=V0[22];  v1100=V[4];  v1010=V[5];
  V0[22]=interpolate_Compare_Middle<ValueType>(2, 1, 1, V[4], V[5], Axis, R, V);

//  0211 0110 0101  v0211=V0[31];  v0110=V[7];  v0101=V[8];
  V0[31]=interpolate_Compare_Middle<ValueType>(0, 2, 1, V[7], V[8], Axis, R, V);

//  1021 0011 1010  v1021=V0[29];  v0011=V[9];  v1010=V[5];
  V0[29]=interpolate_Compare_Middle<ValueType>(1, 0, 2, V[9], V[5], Axis, R, V);

//  1102 1001 0101  v1102=V0[27];  v1001=V[6];  v0101=V[8];
  V0[27]=interpolate_Compare_Middle<ValueType>(1, 1, 0, V[6], V[8], Axis, R, V);


//  1210 0110 1100  v1210=V0[23];  v0110=V[7];  v1100=V[4];
  V0[23]=interpolate_Compare_Middle<ValueType>(1, 2, 1, V[7], V[4], Axis, R, V);

//  0121 0011 0110  v0121=V0[32];  v0011=V[9];  v0110=V[7]; 
  V0[32]=interpolate_Compare_Middle<ValueType>(0, 1, 2, V[9], V[7], Axis, R, V);

//  1012 1001 0011  v1012=V0[30];  v1001=V[6];  v0011=V[9]; 
  V0[30]=interpolate_Compare_Middle<ValueType>(1, 0, 1, V[6], V[9], Axis, R, V);

//  2101 1100 1001  v2101=V0[25];  v1100=V[4];  v1001=V[6];
  V0[25]=interpolate_Compare_Middle<ValueType>(2, 1, 0, V[4], V[6], Axis, R, V);


//  1120 1010 0110  v1120=V0[24];  v1010=V[5];  v0110=V[7];
  V0[24]=interpolate_Compare_Middle<ValueType>(1, 1, 2, V[5], V[7], Axis, R, V);

//  0112 0101 0011  v0112=V0[33];  v0101=V[8];  v0011=V[9]; 
  V0[33]=interpolate_Compare_Middle<ValueType>(0, 1, 1, V[8], V[9], Axis, R, V);

//  2011 1010 1001  v2011=V0[28];  v1010=V[5];  v1001=V[6];
  V0[28]=interpolate_Compare_Middle<ValueType>(2, 0, 1, V[5], V[6], Axis, R, V);

//  1201 0101 1100  v1201=V0[26];  v0101=V[8];  v1100=V[4];
  V0[26]=interpolate_Compare_Middle<ValueType>(1, 2, 0, V[8], V[4], Axis, R, V);


//  1111 v1111=V0[34]; 

//  0) 0110 1001
  if (Axis=='0')
  {
  //  1111 0110 1001  v1111=V0[34];  v0110=V[7];  v1001=V[6];
    V0[34]=interpolate_Compare_Middle<ValueType>(1, 1, 1, V[7], V[6], Axis, R, V);
  }

//  1) 1010 0101
  else if (Axis=='1')
  {
  //  1111 1010 0101  v1111=V0[34];  v1010=V[5];  v0101=V[8];
    V0[34]=interpolate_Compare_Middle<ValueType>(1, 1, 1, V[5], V[8], Axis, R, V);
  }

//  2) 1100 0011
  else
  {
  //  1111 1100 0011  v1111=V0[34];  v1100=V[4];  v0011=V[9];
    V0[34]=interpolate_Compare_Middle<ValueType>(1, 1, 1, V[4], V[9], Axis, R, V);
  };

  return(Interpolate_3D_interpolation_2_4<ValueType>('P', Axis, R, V0));
}

/////////////////////////////////interpolation_2_4///////////////////////////////////////////////////////////

template <typename ValueType>
ValueType Interpolate_3D_interpolation_2_4(
const char How, 
const char Axis,
const ValueType * R,
const ValueType * V)
{
  ValueType x, y, z;
  ValueType xA, yA, zA;
  ValueType xB, yB, zB;
  ValueType xC, yC, zC;
  ValueType xD, yD, zD;
  ValueType v2000, v0200, v0020, v0002;
  ValueType v1100, v0110, v0011, v1001, v1010, v0101;

  ValueType v4000, v0400, v0040, v0004;
  ValueType v3100, v0310, v0031, v1003; 
  ValueType v1300, v0130, v0013, v3001; 
  ValueType v3010, v0301, v1030, v0103;
  ValueType v2200, v0220, v0022, v2002, v2020, v0202;
  ValueType v2110, v0211, v1021, v1102; 
  ValueType v1210, v0121, v1012, v2101; 
  ValueType v1120, v0112, v2011, v1201;
  ValueType v1111;

  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;
  ValueType R0 [15];
  ValueType V0 [10];

  x=R[0];
  y=R[1];
  z=R[2];

  v4000=V[0]; 
  v0400=V[1]; 
  v0040=V[2]; 
  v0004=V[3];

  v3100=V[4];
  v2200=V[5];
  v1300=V[6];

  v3010=V[7];
  v2020=V[8];
  v1030=V[9];

  v3001=V[10];
  v2002=V[11];
  v1003=V[12];
 
  v0310=V[13]; 
  v0220=V[14];
  v0130=V[15];

  v0301=V[16];
  v0202=V[17];
  v0103=V[18];

  v0031=V[19]; 
  v0022=V[20]; 
  v0013=V[21]; 

  v2110=V[22];
  v1210=V[23];
  v1120=V[24];
 
  v2101=V[25];
  v1201=V[26];
  v1102=V[27];

  v2011=V[28];
  v1021=V[29];
  v1012=V[30];

  v0211=V[31]; 
  v0121=V[32]; 
  v0112=V[33]; 

  v1111=V[34];

  R0[0]=R[0];
  R0[1]=R[1];
  R0[2]=R[2];

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     R[3], R[4], R[5],  
                                     R[6], R[7], R[8],  
                                     R[9], R[10], R[11],  
                                     R[12], R[13], R[14], 
                                     &sA, &sB, &sC, &sD);
  sA=(ValueType)2.0*sA;
  sB=(ValueType)2.0*sB;
  sC=(ValueType)2.0*sC;
  sD=(ValueType)2.0*sD;  

  if (sA>(ValueType)1) //1 Тетраэдр 4000 2200 2020 2002
  {
    interpolation_ABCD_Coordinate<ValueType>(4, 0, 0, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xB, &yB, &zB, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);
 
    v2000=v4000;
    v0200=v2200;
    v0020=v2020;
    v0002=v2002;

    v1100=v3100;
    v0110=v2110;
    v0011=v2011;
    v1001=v3001;
    v1010=v3010;
    v0101=v2101;
  }
  else if (sB>(ValueType)1) //2 Тетраэдр 2200 0400 0220 0202
  {
    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 4, 0, 0, 4,
                                  &xB, &yB, &zB, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);
 
    v2000=v2200;
    v0200=v0400;
    v0020=v0220;
    v0002=v0202;

    v1100=v1300;
    v0110=v0310;
    v0011=v0211;
    v1001=v1201;
    v1010=v1210;
    v0101=v0301;
  }
  else if (sC>(ValueType)1) //3 Тетраэдр 2020 0220 0040 0022
  {
    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xB, &yB, &zB, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 4, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);
 
    v2000=v2020;
    v0200=v0220;
    v0020=v0040;
    v0002=v0022;

    v1100=v1120;
    v0110=v0130;
    v0011=v0031;
    v1001=v1021;
    v1010=v1030;
    v0101=v0121;
  }
  else if (sD>(ValueType)1) //4 Тетраэдр 2002 0202 0022 0004
  {
    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xA, &yA, &zA, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xC, &yC, &zC, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 0, 4, 4,
                                  &xD, &yD, &zD, 
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11],  
                                  R[12], R[13], R[14]);
 
    v2000=v2002;
    v0200=v0202;
    v0020=v0022;
    v0002=v0004;

    v1100=v1102;
    v0110=v0112;
    v0011=v0013;
    v1001=v1003;
    v1010=v1012;
    v0101=v0103;
  }
  else //5
  {
    if(Axis=='0') //5.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0) 
        {
//5.1.1 Тетраэдр 0220 2002 2200 2020
          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v0220;
          v0200=v2002;
          v0020=v2200;
          v0002=v2020;

          v1100=v1111;
          v0110=v2101;
          v0011=v2110;
          v1001=v1120;
          v1010=v1210;
          v0101=v2011;
        }
        else 
        {
//5.1.3 Тетраэдр 0220 2002 2200 0202
          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v0220;
          v0200=v2002;
          v0020=v2200;
          v0002=v0202;

          v1100=v1111;
          v0110=v2101;
          v0011=v1201;
          v1001=v0211;
          v1010=v1210;
          v0101=v1102;
        };
      }
      else
      {
        if (s2>(ValueType)0) 
        {
//5.1.4 Тетраэдр 0220 2002 0022 2020
          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v0220;
          v0200=v2002;
          v0020=v0022;
          v0002=v2020;

          v1100=v1111;
          v0110=v1012;
          v0011=v1021;
          v1001=v1120;
          v1010=v0121;
          v0101=v2011;
        }
        else 
        {
//5.1.2 Тетраэдр 0220 2002 0022 0202
          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v0220;
          v0200=v2002;
          v0020=v0022;
          v0002=v0202;

          v1100=v1111;
          v0110=v1012;
          v0011=v0112;
          v1001=v0211;
          v1010=v0121;
          v0101=v1102;
        };
      };
    }
    else if(Axis=='1') //5.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0) 
        {
//5.2.1 Тетраэдр 2020 0202 2200 0220
          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2020;
          v0200=v0202;
          v0020=v2200;
          v0002=v0220;

          v1100=v1111;
          v0110=v1201;
          v0011=v1210;
          v1001=v1120;
          v1010=v2110;
          v0101=v0211;
        }
        else 
        {
//5.2.3 Тетраэдр 2020 0202 2200 2002
          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2020;
          v0200=v0202;
          v0020=v2200;
          v0002=v2002;

          v1100=v1111;
          v0110=v1201;
          v0011=v2101;
          v1001=v2011;
          v1010=v2110;
          v0101=v1102;
        };
      }
      else
      {
        if (s2>(ValueType)0) 
        {
//5.2.2 Тетраэдр 2020 0202 0022 0220
          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2020;
          v0200=v0202;
          v0020=v0022;
          v0002=v0220;

          v1100=v1111;
          v0110=v0112;
          v0011=v0121;
          v1001=v1120;
          v1010=v1021;
          v0101=v0211;
        }
        else 
        {
//5.2.4 Тетраэдр 2020 0202 0022 2002
          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2020;
          v0200=v0202;
          v0020=v0022;
          v0002=v2002;

          v1100=v1111;
          v0110=v0112;
          v0011=v1012;
          v1001=v2011;
          v1010=v1021;
          v0101=v1102;
        };
      };
    }
    else //5.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0) 
        {
//5.3.1 Тетраэдр 2200 0022 0220 2020
          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2200;
          v0200=v0022;
          v0020=v0220;
          v0002=v2020;

          v1100=v1111;
          v0110=v0121;
          v0011=v1120;
          v1001=v2110;
          v1010=v1210;
          v0101=v1021;
        }
        else 
        {
//5.3.2 Тетраэдр 2200 0022 0220 0202
          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2200;
          v0200=v0022;
          v0020=v0220;
          v0002=v0202;

          v1100=v1111;
          v0110=v0121;
          v0011=v0211;
          v1001=v1201;
          v1010=v1210;
          v0101=v0112;
        };
      }
      else
      {
        if (s2>(ValueType)0) 
        {
//5.3.4 Тетраэдр 2200 0022 2002 2020
          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2200;
          v0200=v0022;
          v0020=v2002;
          v0002=v2020;

          v1100=v1111;
          v0110=v1012;
          v0011=v2011;
          v1001=v2110;
          v1010=v2101;
          v0101=v1021;
        }
        else 
        {
//5.3.3 Тетраэдр 2200 0022 2002 0202
          interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                        &xA, &yA, &zA, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                        &xB, &yB, &zB, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                        &xC, &yC, &zC, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);

          interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                        &xD, &yD, &zD, 
                                        R[3], R[4], R[5],  
                                        R[6], R[7], R[8],  
                                        R[9], R[10], R[11],  
                                        R[12], R[13], R[14]);
       
          v2000=v2200;
          v0200=v0022;
          v0020=v2002;
          v0002=v0202;

          v1100=v1111;
          v0110=v1012;
          v0011=v1102;
          v1001=v1201;
          v1010=v2101;
          v0101=v0112;
        };
      };
    };        
  };

  R0[3]=xA;
  R0[4]=yA;
  R0[5]=zA;

  R0[6]=xB;
  R0[7]=yB;
  R0[8]=zB;

  R0[9]=xC;
  R0[10]=yC;
  R0[11]=zC;

  R0[12]=xD;
  R0[13]=yD;
  R0[14]=zD;

  V0[0]=v2000;
  V0[1]=v0200;
  V0[2]=v0020;
  V0[3]=v0002;

  V0[4]=v1100;
  V0[5]=v1010;
  V0[6]=v1001;
  V0[7]=v0110;
  V0[8]=v0101;
  V0[9]=v0011;

  result=Interpolate_3D<ValueType>(2, How, '0', R0, V0);

  return(result);
}

template <typename ValueType>
ValueType interpolation_MinAltitude_3D_interpolation_2_4_axis_1(
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD)
{
  ValueType xAB, yAB, zAB;
  ValueType xBC, yBC, zBC;
  ValueType xCD, yCD, zCD;
  ValueType xDA, yDA, zDA;
  ValueType xAC, yAC, zAC;
  ValueType xBD, yBD, zBD;

  ValueType V;
  ValueType S_A, S_B, S_C, S_D;
  ValueType S_1, S_2;
  ValueType S1, S2, S3;
  ValueType Result;

  V=interpolationCompositionalProduct3D<ValueType>(xC-xD, yC-yD, zC-zD, xA-xD, yA-yD, zA-zD, xB-xD, yB-yD, zB-zD);
  if (V<(ValueType)0)
  {
    V=(ValueType)0-V;
  };

  S_A=interpolationAreaOfTriangle_3D<ValueType>(xB, yB, zB,
                                     xC, yC, zC,
                                     xD, yD, zD);

  S_B=interpolationAreaOfTriangle_3D<ValueType>(xC, yC, zC,
                                     xD, yD, zD,
                                     xA, yA, zA);

  S_C=interpolationAreaOfTriangle_3D<ValueType>(xD, yD, zD,
                                     xA, yA, zA,
                                     xB, yB, zB);

  S_D=interpolationAreaOfTriangle_3D<ValueType>(xA, yA, zA,
                                     xB, yB, zB,
                                     xC, yC, zC);


  S_A=(ValueType)0.25*S_A;
  S_B=(ValueType)0.25*S_B;
  S_C=(ValueType)0.25*S_C;
  S_D=(ValueType)0.25*S_D;

  xAB=(xA+xB)*(ValueType)0.5;
  yAB=(yA+yB)*(ValueType)0.5;
  zAB=(zA+zB)*(ValueType)0.5;

  xBC=(xB+xC)*(ValueType)0.5;
  yBC=(yB+yC)*(ValueType)0.5;
  zBC=(zB+zC)*(ValueType)0.5;

  xCD=(xC+xD)*(ValueType)0.5;
  yCD=(yC+yD)*(ValueType)0.5;
  zCD=(zC+zD)*(ValueType)0.5;

  xDA=(xD+xA)*(ValueType)0.5;
  yDA=(yD+yA)*(ValueType)0.5;
  zDA=(zD+zA)*(ValueType)0.5;

  xAC=(xA+xC)*(ValueType)0.5;
  yAC=(yA+yC)*(ValueType)0.5;
  zAC=(zA+zC)*(ValueType)0.5;

  xBD=(xB+xD)*(ValueType)0.5;
  yBD=(yB+yD)*(ValueType)0.5;
  zBD=(zB+zD)*(ValueType)0.5;

//Axis 1
    S_1=interpolationAreaOfTriangle_3D<ValueType>(xBC, yBC, zBC,
                                       xDA, yDA, zDA,
                                       xAC, yAC, zAC);  

    S_2=interpolationAreaOfTriangle_3D<ValueType>(xBC, yBC, zBC,
                                       xDA, yDA, zDA,
                                       xAB, yAB, zAB);  
    S1=interpolation_Max_2<ValueType>(S_1, S_2);

//Axis 2
//    S_1=interpolationAreaOfTriangle_3D<ValueType>(xAC, yAC, zAC,
//                                       xBD, yBD, zBD,
//                                      xBC, yBC, zBC);  
//
//    S_2=interpolationAreaOfTriangle_3D<ValueType>(xAC, yAC, zAC,
//                                       xBD, yBD, zBD,
//                                       xAB, yAB, zAB);  
//    S2=interpolation_Max_2<ValueType>(S_1, S_2);

//Axis 3
//    S_1=interpolationAreaOfTriangle_3D<ValueType>(xAB, yAB, zAB,
//                                       xCD, yCD, zCD,
//                                       xAC, yAC, zAC);
//  
//    S_2=interpolationAreaOfTriangle_3D(xAB, yAB, zAB,
//                                       xCD, yCD, zCD,
//                                       xBC, yBC, zBC);  
//    S3=interpolation_Max_2<ValueType>(S_1, S_2);

  Result=interpolation_Max_5<ValueType>(S1, S_A, S_B, S_C, S_D);
  Result=V/((ValueType)8.0*Result);

  return(Result);
}


template <typename ValueType>
void interpolation_ABCD_Coordinate(const int a, const int b, const int c, const int d, const int N,
                                   ValueType * X, ValueType * Y, ValueType * Z, 
                                   const ValueType xA, const ValueType yA, const ValueType zA,  
                                   const ValueType xB, const ValueType yB, const ValueType zB,  
                                   const ValueType xC, const ValueType yC, const ValueType zC,  
                                   const ValueType xD, const ValueType yD, const ValueType zD)

{
  (*X)=xD+((((ValueType)a)/((ValueType)N))*(xA-xD)+(((ValueType)b)/((ValueType)N))*(xB-xD)+(((ValueType)c)/((ValueType)N))*(xC-xD));
  (*Y)=yD+((((ValueType)a)/((ValueType)N))*(yA-yD)+(((ValueType)b)/((ValueType)N))*(yB-yD)+(((ValueType)c)/((ValueType)N))*(yC-yD));
  (*Z)=zD+((((ValueType)a)/((ValueType)N))*(zA-zD)+(((ValueType)b)/((ValueType)N))*(zB-zD)+(((ValueType)c)/((ValueType)N))*(zC-zD));

  return;
}



template <typename ValueType>
ValueType Interpolation_MinAltitude_3D_Swift_interpolation_2_4(char * Axis, const ValueType * R)
{
  ValueType Result, Result0;
  ValueType R_A, R_B, R_C, R_D, R_1, R_2, R_3, R_4;
  int a, b, c, d;
  ValueType xA, yA, zA;  
  ValueType xB, yB, zB;  
  ValueType xC, yC, zC;  
  ValueType xD, yD, zD;



// Вычисляем по сути оптимальную ось для тетраэдра в целом
  Result=interpolation_MinAltitude_3D<ValueType>(Axis, 4,
                                      R[0], R[ 1], R[ 2],
                                      R[3], R[ 4], R[ 5],
                                      R[6], R[ 7], R[ 8],
                                      R[9], R[10], R[11]); 

// Вычисляем минимальную высоту для каждого из 8 маленьких тетраэдров с учетом извесной оси. В каждом из них используется 1ая ось. 

//1. Тетраэдр A 4000 2200 2020 2002

  interpolation_ABCD_Coordinate<ValueType>(4, 0, 0, 0, 4,
                                &xA, &yA, &zA, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                &xB, &yB, &zB, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                &xC, &yC, &zC, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                &xD, &yD, &zD, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  R_A=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                            xB, yB, zB,
                                                            xC, yC, zC,
                                                            xD, yD, zD);


//2. Тетраэдр B 2200 0400 0220 0202

  interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                &xA, &yA, &zA, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 4, 0, 0, 4,
                                &xB, &yB, &zB, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                &xC, &yC, &zC, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                &xD, &yD, &zD, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  R_B=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                            xB, yB, zB,
                                                            xC, yC, zC,
                                                            xD, yD, zD);

//3. Тетраэдр C 2020 0220 0040 0022

  interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                &xA, &yA, &zA, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                &xB, &yB, &zB, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 0, 4, 0, 4,
                                &xC, &yC, &zC, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                &xD, &yD, &zD, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  R_C=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                            xB, yB, zB,
                                                            xC, yC, zC,
                                                            xD, yD, zD);

//4. Тетраэдр D 2002 0202 0022 0004

  interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                &xA, &yA, &zA, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                &xB, &yB, &zB, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                &xC, &yC, &zC, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  interpolation_ABCD_Coordinate<ValueType>(0, 0, 0, 4, 4,
                                &xD, &yD, &zD, 
                                R[0], R[1], R[2],  
                                R[3], R[4], R[5],  
                                R[6], R[7], R[8],  
                                R[9], R[10], R[11]);

  R_D=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                            xB, yB, zB,
                                                            xC, yC, zC,
                                                            xD, yD, zD);



// Ось 5.1.
  if ((*Axis)=='0')
  {
  //5.1.1. Тетраэдр 0220 2002 2200 2020

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_1=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.1.2. Тетраэдр 0220 2002 0022 0202

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_2=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.1.3. Тетраэдр 0220 2002 2200 0202

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_3=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.1.4. Тетраэдр 0220 2002 0022 2020

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_4=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);
  }
// Ось 5.2.
  else if ((*Axis)=='1')
  {
  //5.2.1. Тетраэдр 2020 0202 2200 0220

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_1=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.2.2. Тетраэдр 2020 0202 0022 0220

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_2=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.2.3. Тетраэдр 2020 0202 2200 2002

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_3=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.2.4. Тетраэдр 2020 0202 0022 2002

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_4=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);
  }

// Ось 5.3.
  else if ((*Axis)=='2')
  {
  //5.3.1. Тетраэдр 2200 0022 0220 2020

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_1=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.3.2. Тетраэдр 2200 0022 0220 0202

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 2, 0, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_2=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.3.3. Тетраэдр 2200 0022 2002 0202

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 2, 0, 2, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_3=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);

  //5.3.4. Тетраэдр 2200 0022 2002 2020

    interpolation_ABCD_Coordinate<ValueType>(2, 2, 0, 0, 4,
                                  &xA, &yA, &zA, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(0, 0, 2, 2, 4,
                                  &xB, &yB, &zB, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 0, 2, 4,
                                  &xC, &yC, &zC, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    interpolation_ABCD_Coordinate<ValueType>(2, 0, 2, 0, 4,
                                  &xD, &yD, &zD, 
                                  R[0], R[1], R[2],  
                                  R[3], R[4], R[5],  
                                  R[6], R[7], R[8],  
                                  R[9], R[10], R[11]);

    R_4=interpolation_MinAltitude_3D_interpolation_2_4_axis_1<ValueType>(xA, yA, zA,
                                                              xB, yB, zB,
                                                              xC, yC, zC,
                                                              xD, yD, zD);
  };

// Далее функцией interpolation_Min_4 сравниваем пополовинно. Далее сравниваем 2 результата и меньший из них уходит в ответ.
  Result=interpolation_Min_4<ValueType>(R_A, R_B, R_C, R_D);
  Result0=interpolation_Min_4<ValueType>(R_1, R_2, R_3, R_4);
  if (Result0<Result)
  {
    Result=Result0;
  };
  return(Result);
}



///////////////////Перестановки для контактов//////////////////

int GetRelativeOrientationSide(int * firstNodes, int * secondNodes)
{
  int i, j, k;
  int a, b, c;
  int R;

  if (firstNodes[0]<firstNodes[1])
  {
    if (firstNodes[1]<firstNodes[2])
    {
      //012
      i=0;
      j=1;
      k=2;
    }
    else
    {
      if (firstNodes[0]<firstNodes[2])
      {
        //021
        i=0;
        j=2;
        k=1;
      }
      else
      {
        //120
        i=2;
        j=0;
        k=1;
      };
    };
  }
  else
  {
    if (firstNodes[1]>firstNodes[2])
    {
      //210
      i=2;
      j=1;
      k=0;
    }
    else
    {
      if (firstNodes[0]<firstNodes[2])
      {
        //102
        i=1;
        j=0;
        k=2;
      }
      else
      {
        //201
        i=1;
        j=2;
        k=0;
      };
    };
  };  

  a=secondNodes[i];
  b=secondNodes[j];
  c=secondNodes[k];

  if (a<b)
  {
    if (b<c)
    {
      R=12;
    }
    else
    {
      if (a<c)
      {
        R=1021;
      }
      else
      {
        R=120;
      };
    };
  }
  else
  {
    if (b>c)
    {
      R=1210;
    }
    else
    {
      if (a<c)
      {
        R=1102;
      }
      else
      {
        R=201;
      };
    };
  };  

  return(R);
}

int GetCorrespondingIndexSide(int relativeOrientation, int firstIndex, int number)
{
  int L;
  int P[3];
  int N;
  int R;

  N=relativeOrientation;
  L=N/1000;

  N=N-1000*L;
  P[0]=N/100;

  N=N-100*P[0];
  P[1]=N/10;

  N=N-10*P[1];
  P[2]=N;

  if (number<4)
  {
    R=firstIndex;
  }

  else if (number==4)
  {
    R=P[firstIndex];
  }

  else if (number==5)
  {
    if (firstIndex<3)
    {
      R=P[firstIndex];
    }
    else
    {
      if (L==0)
      {
        R=P[firstIndex-3]+3;
      }
      else
      {
        if (firstIndex==3)
        {
          R=P[1]+3;
        }
        else if (firstIndex==4)
        {
          R=P[2]+3;
        }
        else if (firstIndex==5)
        {
          R=P[0]+3;
        };
      };
    };
  };
  return(R);
}

int GetRelativeOrientationEdge(int * firstNodes, int * secondNodes)
{
  int R;
  if (firstNodes[0]<firstNodes[1])
  {
    if (secondNodes[0]<secondNodes[1])
    {
      R=0;
    }
    else
    {
      R=1;
    };
  }
  else
  {
    if (secondNodes[0]<secondNodes[1])
    {
      R=1;
    }
    else
    {
      R=0;
    };
  };
  return(R);
}

int GetCorrespondingIndexEdge(int relativeOrientation, int firstIndex, int number)
{
  int R;
  if ((number>2)&&(relativeOrientation==1))
  {
    R=number-2-firstIndex;
  }
  else 
  {
    R=firstIndex;
  };
  return(R);
}

///////////////////////MaxAltitude////////////////////////////////

template <typename ValueType>
ValueType Interpolation_MaxAltitude_3D_Swift(char Axis, const int N, const ValueType * R)
{
  ValueType xA, yA, zA;
  ValueType xB, yB, zB;
  ValueType xC, yC, zC;
  ValueType xD, yD, zD;

  ValueType xAB, yAB, zAB;
  ValueType xBC, yBC, zBC;
  ValueType xCD, yCD, zCD;
  ValueType xDA, yDA, zDA;
  ValueType xAC, yAC, zAC;
  ValueType xBD, yBD, zBD;

  ValueType V;
  ValueType S_A, S_B, S_C, S_D;
  ValueType S_1, S_2;
  ValueType Result;

  xA=R[0];
  yA=R[1];
  zA=R[2];
  xB=R[3];
  yB=R[4];
  zB=R[5];
  xC=R[6];
  yC=R[7];
  zC=R[8];
  xD=R[9];
  yD=R[10];
  zD=R[11];

  V=interpolationCompositionalProduct3D<ValueType>(xC-xD, yC-yD, zC-zD, xA-xD, yA-yD, zA-zD, xB-xD, yB-yD, zB-zD);
  if (V<(ValueType)0)
  {
    V=(ValueType)0-V;
  };

  S_A=interpolationAreaOfTriangle_3D<ValueType>(xB, yB, zB,
                                     xC, yC, zC,
                                     xD, yD, zD);

  S_B=interpolationAreaOfTriangle_3D<ValueType>(xC, yC, zC,
                                     xD, yD, zD,
                                     xA, yA, zA);

  S_C=interpolationAreaOfTriangle_3D<ValueType>(xD, yD, zD,
                                     xA, yA, zA,
                                     xB, yB, zB);

  S_D=interpolationAreaOfTriangle_3D<ValueType>(xA, yA, zA,
                                     xB, yB, zB,
                                     xC, yC, zC);

  if (N==(int)1)
  {
    Result=interpolation_Min_4<ValueType>(S_A, S_B, S_C, S_D);
    Result=V/Result;
  }
  else
  {
    S_A=(ValueType)0.25*S_A;
    S_B=(ValueType)0.25*S_B;
    S_C=(ValueType)0.25*S_C;
    S_D=(ValueType)0.25*S_D;

    xAB=(xA+xB)*(ValueType)0.5;
    yAB=(yA+yB)*(ValueType)0.5;
    zAB=(zA+zB)*(ValueType)0.5;

    xBC=(xB+xC)*(ValueType)0.5;
    yBC=(yB+yC)*(ValueType)0.5;
    zBC=(zB+zC)*(ValueType)0.5;

    xCD=(xC+xD)*(ValueType)0.5;
    yCD=(yC+yD)*(ValueType)0.5;
    zCD=(zC+zD)*(ValueType)0.5;

    xDA=(xD+xA)*(ValueType)0.5;
    yDA=(yD+yA)*(ValueType)0.5;
    zDA=(zD+zA)*(ValueType)0.5;

    xAC=(xA+xC)*(ValueType)0.5;
    yAC=(yA+yC)*(ValueType)0.5;
    zAC=(zA+zC)*(ValueType)0.5;

    xBD=(xB+xD)*(ValueType)0.5;
    yBD=(yB+yD)*(ValueType)0.5;
    zBD=(zB+zD)*(ValueType)0.5;

    if (Axis=='0')
    {
      S_1=interpolationAreaOfTriangle_3D<ValueType>(xBC, yBC, zBC,
                                         xDA, yDA, zDA,
                                         xAC, yAC, zAC);  
  
      S_2=interpolationAreaOfTriangle_3D<ValueType>(xBC, yBC, zBC,
                                         xDA, yDA, zDA,
                                         xAB, yAB, zAB);  
    }
    else if (Axis=='1')
    {
      S_1=interpolationAreaOfTriangle_3D<ValueType>(xAC, yAC, zAC,
                                         xBD, yBD, zBD,
                                         xBC, yBC, zBC);  

      S_2=interpolationAreaOfTriangle_3D<ValueType>(xAC, yAC, zAC,
                                         xBD, yBD, zBD,
                                         xAB, yAB, zAB);  
    }
    else
    {
      S_1=interpolationAreaOfTriangle_3D<ValueType>(xAB, yAB, zAB,
                                         xCD, yCD, zCD,
                                         xAC, yAC, zAC);
  
      S_2=interpolationAreaOfTriangle_3D<ValueType>(xAB, yAB, zAB,
                                         xCD, yCD, zCD,
                                         xBC, yBC, zBC);  
    };
    Result=interpolation_Min_6<ValueType>(S_1, S_2, S_A, S_B, S_C, S_D);
    Result=V/((ValueType)4.0*(ValueType)N*Result);
  };
  return(Result);
}

template <typename ValueType>
ValueType interpolation_Min_4(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4)
{
  ValueType R, R1;

  if (s1>s2)
  {
    R=s2;
  }
  else
  {
    R=s1;
  };

  if (s3>s4)
  {
    R1=s4;
  }
  else
  {
    R1=s3;
  };

  if (R>R1)
  {
    R=R1;
  };

  return(R);
}


template <typename ValueType>
ValueType interpolation_Min_6(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4, const ValueType s5, const ValueType s6)
{
  ValueType R, R1, R2;

  if (s1>s2)
  {
    R=s2;
  }
  else
  {
    R=s1;
  };

  if (s3>s4)
  {
    R1=s4;
  }
  else
  {
    R1=s3;
  };

  if (s5>s6)
  {
    R2=s6;
  }
  else
  {
    R2=s5;
  };

  if (R>R1)
  {
    R=R1;
  };

  if (R>R2)
  {
    R=R2;
  };

  return(R);
}

///////////////////////Permitaions////////////////////////////////////////////

template <typename ValueType>
void InterpolationPermitation_Edge(const int N, const int * P, ValueType * V)
//N>1, P[2], V[N-1]
{
  if (N>2)
  {
    if(P[0]>P[1])
    {
      int i, J;
      ValueType V1;

      J=(N-1)/2;
      for (i=0; (i<J); i++)
      {
        V1=V[N-2-i];
        V[N-2-i]=V[i];
        V[i]=V1;
      };
    };
  };
  return;
}

template <typename ValueType>
void InterpolationPermitation_Side(const int N, const int * P, ValueType * V)
//N>2, P[3], V[(N-1)*(N-2)/2]
{
  if (N>3)
  {
    int I, j, k;
    int P1[2];
    int I1[2];
    ValueType V1[6];

    if (P[0]<P[1])
    {
      I=0;
      P1[0]=P[0];
    }
    else
    {
      I=1;
      P1[0]=P[1];
    };
    if (P[2]<P1[0])
    {
      I=2;
    };
    
    for(j=0, k=0; (j<3); j++)
    {
      if (j!=I)
      {
        P1[k]=P[j];
        I1[k]=j;
        k++;
      };
    };

    V1[0]=V[I];
    if (N==5)
    {
      if (I==2)
      {
        V1[4]=V[3];
      }
      else
      {
        V1[4]=V[4+I];
      };
    };

    if (P1[0]<P1[1])
    {
      I=I1[0];
      j=I1[1];
    }
    else
    {
      I=I1[1];
      j=I1[0];
    };
    
    V1[1]=V[I];
    if (N==5)
    {
      if (I==2)
      {
        V1[5]=V[3];
      }
      else
      {
        V1[5]=V[4+I];
      };
    };
    
    V1[2]=V[j];
    if (N==5)
    {
      if (j==2)
      {
        V1[3]=V[3];
      }
      else
      {
        V1[3]=V[4+j];
      };
    };

    for(k=0; (k<((N-1)*(N-2)/2)); k++)
    {
      V[k]=V1[k];
    };
  };
  return;
}

template <typename ValueType>
void InterpolationPermitation_Volume(const int N, const int * P, ValueType * V)
//N>3, P[4], V[(N-1)*(N-2)*(N-3)/6]
{
  if (N==5)
  {
    int I, j, k;
    int P1[3];
    int I1[3];
    ValueType V1[4];

    if (P[0]<P[1])
    {
      I=0;
      P1[0]=P[0];
    }
    else
    {
      I=1;
      P1[0]=P[1];
    };
    if (P[2]<P[3])
    {
      j=2;
      P1[1]=P[2];
    }
    else
    {
      j=3;
      P1[1]=P[3];
    };
    if (P1[1]<P1[0])
    {
      I=j;
    };

    for(j=0, k=0; (j<4); j++)
    {
      if (j!=I)
      {
        P1[k]=P[j];
        I1[k]=j;
        k++;
      };
    };
    V1[0]=V[I];

    if (P1[0]<P1[1])
    {
      I=I1[0];
      j=P1[0];
    }
    else
    {
      I=I1[1];
      j=P1[1];
    };
    if (P1[2]<j)
    {
      I=I1[2];
    };

    for(j=0, k=0; (j<3); j++)
    {
      if (j!=I)
      {
        P1[k]=P1[j];
        I1[k]=I1[j];
        k++;
      };
    };
    V1[1]=V[I];

    if (P1[0]<P1[1])
    {
      I=I1[0];
      j=I1[1];
    }
    else
    {
      I=I1[1];
      j=I1[0];
    };
    
    V1[2]=V[I];
    V1[3]=V[j];

    for(k=0; (k<4); k++)
    {
      V[k]=V1[k];
    };
  };
  return;
}

///////////////////////Courant////////////////////////////////////////////

template <typename ValueType>
ValueType interpolationAreaOfTriangle_3D(
ValueType xA, ValueType yA, ValueType zA,
ValueType xB, ValueType yB, ValueType zB,
ValueType xC, ValueType yC, ValueType zC)
{
  ValueType Result;

  xB=xB-xA;
  yB=yB-yA;
  zB=zB-zA;

  xC=xC-xA;
  yC=yC-yA;
  zC=zC-zA;

  xA=yB*zC-zB*yC;
  yA=zB*xC-xB*zC;
  zA=xB*yC-yB*xC;

  Result=xA*xA+yA*yA+zA*zA;
  Result=(ValueType)sqrt((double)Result);
  return(Result);
}

template <typename ValueType>
ValueType interpolation_Max_2(const ValueType s1, const ValueType s2)
{
  if (s1<s2)
  {
    return(s2);
  };
  return(s1);
}

template <typename ValueType>
ValueType interpolation_MinAndNumber_3(char * a,
const ValueType s0, const ValueType s1, const ValueType s2)
{
  ValueType Result;
  if (s0<s1)
  {
    Result=s0;
    (*a)='0';
  }
  else
  {
    Result=s1;
    (*a)='1';
  };
  if (s2<Result)
  {
    Result=s2;
    (*a)='2';
  };
  return(Result);  
}

template <typename ValueType>
ValueType interpolation_Max_4(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4)
{
  ValueType R, R1;

  if (s1<s2)
  {
    R=s2;
  }
  else
  {
    R=s1;
  };

  if (s3<s4)
  {
    R1=s4;
  }
  else
  {
    R1=s3;
  };

  if (R<R1)
  {
    R=R1;
  };

  return(R);
}

template <typename ValueType>
ValueType interpolation_Max_5(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4, const ValueType s5)
{
  ValueType R, R1;

  if (s1<s2)
  {
    R=s2;
  }
  else
  {
    R=s1;
  };

  if (s4<s5)
  {
    R1=s5;
  }
  else
  {
    R1=s4;
  };

  if (R<s3)
  {
    R=s3;
  };

  if (R<R1)
  {
    R=R1;
  };

  return(R);
}

template <typename ValueType>
ValueType interpolation_MinAltitude_3D(char * Axis, const int N,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD)
{
  ValueType xAB, yAB, zAB;
  ValueType xBC, yBC, zBC;
  ValueType xCD, yCD, zCD;
  ValueType xDA, yDA, zDA;
  ValueType xAC, yAC, zAC;
  ValueType xBD, yBD, zBD;

  ValueType V;
  ValueType S_A, S_B, S_C, S_D;
  ValueType S_1, S_2;
  ValueType S1, S2, S3;
  ValueType Result;

  V=interpolationCompositionalProduct3D<ValueType>(xC-xD, yC-yD, zC-zD, xA-xD, yA-yD, zA-zD, xB-xD, yB-yD, zB-zD);
  if (V<(ValueType)0)
  {
    V=(ValueType)0-V;
  };

  S_A=interpolationAreaOfTriangle_3D<ValueType>(xB, yB, zB,
                                     xC, yC, zC,
                                     xD, yD, zD);

  S_B=interpolationAreaOfTriangle_3D<ValueType>(xC, yC, zC,
                                     xD, yD, zD,
                                     xA, yA, zA);

  S_C=interpolationAreaOfTriangle_3D<ValueType>(xD, yD, zD,
                                     xA, yA, zA,
                                     xB, yB, zB);

  S_D=interpolationAreaOfTriangle_3D<ValueType>(xA, yA, zA,
                                     xB, yB, zB,
                                     xC, yC, zC);

  if (N==(int)1)
  {
    (*Axis)='0';
    Result=interpolation_Max_4<ValueType>(S_A, S_B, S_C, S_D);
    Result=V/Result;
  }
  else
  {
    S_A=(ValueType)0.25*S_A;
    S_B=(ValueType)0.25*S_B;
    S_C=(ValueType)0.25*S_C;
    S_D=(ValueType)0.25*S_D;

    xAB=(xA+xB)*(ValueType)0.5;
    yAB=(yA+yB)*(ValueType)0.5;
    zAB=(zA+zB)*(ValueType)0.5;

    xBC=(xB+xC)*(ValueType)0.5;
    yBC=(yB+yC)*(ValueType)0.5;
    zBC=(zB+zC)*(ValueType)0.5;

    xCD=(xC+xD)*(ValueType)0.5;
    yCD=(yC+yD)*(ValueType)0.5;
    zCD=(zC+zD)*(ValueType)0.5;

    xDA=(xD+xA)*(ValueType)0.5;
    yDA=(yD+yA)*(ValueType)0.5;
    zDA=(zD+zA)*(ValueType)0.5;

    xAC=(xA+xC)*(ValueType)0.5;
    yAC=(yA+yC)*(ValueType)0.5;
    zAC=(zA+zC)*(ValueType)0.5;

    xBD=(xB+xD)*(ValueType)0.5;
    yBD=(yB+yD)*(ValueType)0.5;
    zBD=(zB+zD)*(ValueType)0.5;

//Axis 1
    S_1=interpolationAreaOfTriangle_3D<ValueType>(xBC, yBC, zBC,
                                       xDA, yDA, zDA,
                                       xAC, yAC, zAC);  

    S_2=interpolationAreaOfTriangle_3D<ValueType>(xBC, yBC, zBC,
                                       xDA, yDA, zDA,
                                       xAB, yAB, zAB);  
    S1=interpolation_Max_2<ValueType>(S_1, S_2);

//Axis 2
    S_1=interpolationAreaOfTriangle_3D<ValueType>(xAC, yAC, zAC,
                                       xBD, yBD, zBD,
                                       xBC, yBC, zBC);  

    S_2=interpolationAreaOfTriangle_3D<ValueType>(xAC, yAC, zAC,
                                       xBD, yBD, zBD,
                                       xAB, yAB, zAB);  
    S2=interpolation_Max_2<ValueType>(S_1, S_2);

//Axis 3
    S_1=interpolationAreaOfTriangle_3D<ValueType>(xAB, yAB, zAB,
                                       xCD, yCD, zCD,
                                       xAC, yAC, zAC);
  
    S_2=interpolationAreaOfTriangle_3D<ValueType>(xAB, yAB, zAB,
                                       xCD, yCD, zCD,
                                       xBC, yBC, zBC);  
    S3=interpolation_Max_2<ValueType>(S_1, S_2);

    Result=interpolation_MinAndNumber_3<ValueType>(Axis, S1, S2, S3);
    Result=interpolation_Max_5<ValueType>(Result, S_A, S_B, S_C, S_D);
    Result=V/((ValueType)4.0*(ValueType)N*Result);
  };
  return(Result);
}

template <typename ValueType>
ValueType Interpolation_MinAltitude_3D_Swift(char * Axis, const int N, const ValueType * R)
{
  return(interpolation_MinAltitude_3D<ValueType>(Axis, N,
                                      R[0], R[ 1], R[ 2],
                                      R[3], R[ 4], R[ 5],
                                      R[6], R[ 7], R[ 8],
                                      R[9], R[10], R[11]));  
}

////////////////общая функция для тетраэдров////////////////

template <typename ValueType>
ValueType Interpolate_3D(
const int N, 
const char How, 
const char Axis,
const ValueType * R,
const ValueType * V)
{
  ValueType Result=(ValueType)0;
  if(N==1)
  {
    Result=Interpolate_3D_Linear<ValueType>(R[ 0], R[ 1], R[ 2],
                                            R[ 3], R[ 4], R[ 5],
                                            R[ 6], R[ 7], R[ 8],
                                            R[ 9], R[10], R[11],
                                            R[12], R[13], R[14],
                                            V[0], V[1], V[2], V[3]);
  }
  else if((N==2)&&(How=='P'))
  {
    Result=Interpolate_3D_Power_2<ValueType>(R[ 0], R[ 1], R[ 2],
                                             R[ 3], R[ 4], R[ 5],
                                             R[ 6], R[ 7], R[ 8],
                                             R[ 9], R[10], R[11],
                                             R[12], R[13], R[14],
                                             V[0], V[1], V[2], V[3],
                                             V[4], V[7], V[9], V[6], V[5], V[8]);
  }
  else if((N==2)&&(How=='L'))
  {
    Result=Interpolate_3D_Linear_2<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                              R[ 3], R[ 4], R[ 5],
                                              R[ 6], R[ 7], R[ 8],
                                              R[ 9], R[10], R[11],
                                              R[12], R[13], R[14],
                                              V[0], V[1], V[2], V[3],
                                              V[4], V[7], V[9], V[6], V[5], V[8]);
  }
  else if((N==2)&&(How=='M'))
  {
    Result=Interpolate_3D_Mono_2<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                            R[ 3], R[ 4], R[ 5],
                                            R[ 6], R[ 7], R[ 8],
                                            R[ 9], R[10], R[11],
                                            R[12], R[13], R[14],
                                            V[0], V[1], V[2], V[3],
                                            V[4], V[7], V[9], V[6], V[5], V[8]);
  }
  else if((N==3)&&(How=='P'))
  {
    Result=Interpolate_3D_Power_3<ValueType>(R[ 0], R[ 1], R[ 2],
                                             R[ 3], R[ 4], R[ 5],
                                             R[ 6], R[ 7], R[ 8],
                                             R[ 9], R[10], R[11],
                                             R[12], R[13], R[14],
                                             V[ 0], V[ 1], V[ 2], V[ 3],
                                             V[ 4], V[10], V[14], V[ 9],
                                             V[ 5], V[11], V[15], V[ 8],
                                             V[ 6], V[12], V[ 7], V[13],
                                             V[16], V[19], V[18], V[17]);
  }
  else if((N==3)&&(How=='L'))
  {
    Result=Interpolate_3D_Linear_3<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                              R[ 3], R[ 4], R[ 5],
                                              R[ 6], R[ 7], R[ 8],
                                              R[ 9], R[10], R[11],
                                              R[12], R[13], R[14],
                                              V[ 0], V[ 1], V[ 2], V[ 3],
                                              V[ 4], V[10], V[14], V[ 9],
                                              V[ 5], V[11], V[15], V[ 8],
                                              V[ 6], V[12], V[ 7], V[13],
                                              V[16], V[19], V[18], V[17]);
  }
  else if((N==3)&&(How=='M'))
  {
    Result=Interpolate_3D_Mono_3<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                            R[ 3], R[ 4], R[ 5],
                                            R[ 6], R[ 7], R[ 8],
                                            R[ 9], R[10], R[11],
                                            R[12], R[13], R[14],
                                            V[ 0], V[ 1], V[ 2], V[ 3],
                                            V[ 4], V[10], V[14], V[ 9],
                                            V[ 5], V[11], V[15], V[ 8],
                                            V[ 6], V[12], V[ 7], V[13],
                                            V[16], V[19], V[18], V[17]);
  }
  else if((N==4)&&(How=='P'))
  {
    Result=Interpolate_3D_Power_4<ValueType>(R[ 0], R[ 1], R[ 2],
                                             R[ 3], R[ 4], R[ 5],
                                             R[ 6], R[ 7], R[ 8],
                                             R[ 9], R[10], R[11],
                                             R[12], R[13], R[14],
                                             V[ 0], V[ 1], V[ 2], V[ 3],
                                             V[ 4], V[13], V[19], V[12],
                                             V[ 6], V[15], V[21], V[10],
                                             V[ 7], V[16], V[ 9], V[18],
                                             V[ 5], V[14], V[20], V[11], V[ 8], V[17],
                                             V[22], V[31], V[29], V[27],
                                             V[23], V[32], V[30], V[25],
                                             V[24], V[33], V[28], V[26],
                                             V[34]);
  }
  else if((N==4)&&(How=='L'))
  {
    Result=Interpolate_3D_Linear_4<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                              R[ 3], R[ 4], R[ 5],
                                              R[ 6], R[ 7], R[ 8],
                                              R[ 9], R[10], R[11],
                                              R[12], R[13], R[14],
                                              V[ 0], V[ 1], V[ 2], V[ 3],
                                              V[ 4], V[13], V[19], V[12],
                                              V[ 6], V[15], V[21], V[10],
                                              V[ 7], V[16], V[ 9], V[18],
                                              V[ 5], V[14], V[20], V[11], V[ 8], V[17],
                                              V[22], V[31], V[29], V[27],
                                              V[23], V[32], V[30], V[25],
                                              V[24], V[33], V[28], V[26],
                                              V[34]);
  }
  else if((N==4)&&(How=='M'))
  {
    Result=Interpolate_3D_Mono_4<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                            R[ 3], R[ 4], R[ 5],
                                            R[ 6], R[ 7], R[ 8],
                                            R[ 9], R[10], R[11],
                                            R[12], R[13], R[14],
                                            V[ 0], V[ 1], V[ 2], V[ 3],
                                            V[ 4], V[13], V[19], V[12],
                                            V[ 6], V[15], V[21], V[10],
                                            V[ 7], V[16], V[ 9], V[18],
                                            V[ 5], V[14], V[20], V[11], V[ 8], V[17],
                                            V[22], V[31], V[29], V[27],
                                            V[23], V[32], V[30], V[25],
                                            V[24], V[33], V[28], V[26],
                                            V[34]);
  }
  else if((N==5)&&(How=='P'))
  {
    Result=Interpolate_3D_Power_5<ValueType>(R[ 0], R[ 1], R[ 2],
                                             R[ 3], R[ 4], R[ 5],
                                             R[ 6], R[ 7], R[ 8],
                                             R[ 9], R[10], R[11],
                                             R[12], R[13], R[14],
                                             V[ 0], V[ 1], V[ 2], V[ 3],
                                             V[ 4], V[16], V[24], V[15],
                                             V[ 7], V[19], V[27], V[12],
                                             V[ 8], V[20], V[11], V[23],
                                             V[ 5], V[17], V[25], V[14],
                                             V[ 6], V[18], V[26], V[13],
                                             V[ 9], V[21], V[10], V[22],
                                             V[28], V[46], V[41], V[36],
                                             V[29], V[47], V[42], V[34],
                                             V[30], V[48], V[40], V[35],
                                             V[31], V[49], V[44], V[39],
                                             V[32], V[50], V[45], V[37],
                                             V[33], V[51], V[43], V[38],
                                             V[52], V[53], V[54], V[55]);
  }
  else if((N==5)&&(How=='L'))
  {
    Result=Interpolate_3D_Linear_5<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                              R[ 3], R[ 4], R[ 5],
                                              R[ 6], R[ 7], R[ 8],
                                              R[ 9], R[10], R[11],
                                              R[12], R[13], R[14],
                                              V[ 0], V[ 1], V[ 2], V[ 3],
                                              V[ 4], V[16], V[24], V[15],
                                              V[ 7], V[19], V[27], V[12],
                                              V[ 8], V[20], V[11], V[23],
                                              V[ 5], V[17], V[25], V[14],
                                              V[ 6], V[18], V[26], V[13],
                                              V[ 9], V[21], V[10], V[22],
                                              V[28], V[46], V[41], V[36],
                                              V[29], V[47], V[42], V[34],
                                              V[30], V[48], V[40], V[35],
                                              V[31], V[49], V[44], V[39],
                                              V[32], V[50], V[45], V[37],
                                              V[33], V[51], V[43], V[38],
                                              V[52], V[53], V[54], V[55]);
  }
  else if((N==5)&&(How=='M'))
  {
    Result=Interpolate_3D_Mono_5<ValueType>(Axis, R[ 0], R[ 1], R[ 2],
                                            R[ 3], R[ 4], R[ 5],
                                            R[ 6], R[ 7], R[ 8],
                                            R[ 9], R[10], R[11],
                                            R[12], R[13], R[14],
                                            V[ 0], V[ 1], V[ 2], V[ 3],
                                            V[ 4], V[16], V[24], V[15],
                                            V[ 7], V[19], V[27], V[12],
                                            V[ 8], V[20], V[11], V[23],
                                            V[ 5], V[17], V[25], V[14],
                                            V[ 6], V[18], V[26], V[13],
                                            V[ 9], V[21], V[10], V[22],
                                            V[28], V[46], V[41], V[36],
                                            V[29], V[47], V[42], V[34],
                                            V[30], V[48], V[40], V[35],
                                            V[31], V[49], V[44], V[39],
                                            V[32], V[50], V[45], V[37],
                                            V[33], V[51], V[43], V[38],
                                            V[52], V[53], V[54], V[55]);
  };
  return(Result);
}

//////////////функции для нахождения радиусов опорных точек//////////////

template <typename ValueType>
void Interpolation_3D_Get_Edge(
ValueType * Result, 
const int N, 
const ValueType * R, 
const int I)
{
  ValueType b;
  b=((ValueType)(I+1))/((ValueType)N);
  Result[0]=R[0]+(R[3]-R[0])*b;
  Result[1]=R[1]+(R[4]-R[1])*b;
  Result[2]=R[2]+(R[5]-R[2])*b;
  return;
}

template <typename ValueType>
void Interpolation_3D_Get_Side(
ValueType * Result, 
const int N, 
const ValueType * R, 
const int I)
{
  ValueType b=(ValueType)0;
  ValueType c=(ValueType)0;
  if(N==3) //0
  {
    b=(ValueType)(1.0/3.0);//1
    c=(ValueType)(1.0/3.0);//1
  }
  else if(N==4) //0,1,2
  {
    if(I==0)
    {
      b=(ValueType)0.25;//1
      c=(ValueType)0.25;//1
    }
    else if(I==1)
    {
      b=(ValueType)0.5; //2
      c=(ValueType)0.25;//1
    }
    else
    {
      b=(ValueType)0.25;//1
      c=(ValueType)0.5; //2
    };
  }
  else if(N==5) //0,1,2,3,4,5
  {
    if(I==0)
    {
      b=(ValueType)0.2;//1
      c=(ValueType)0.2;//1
    }
    else if(I==1)
    {
      b=(ValueType)0.6;//3
      c=(ValueType)0.2;//1
    }
    else if(I==2)
    {
      b=(ValueType)0.2;//1
      c=(ValueType)0.6;//3
    }
    else if(I==3)
    {
      b=(ValueType)0.4;//2
      c=(ValueType)0.2;//1
    }
    else if(I==4)
    {
      b=(ValueType)0.4;//2
      c=(ValueType)0.4;//2
    }
    else
    {
      b=(ValueType)0.2;//1
      c=(ValueType)0.4;//2
    };
  };
  Result[0]=R[0]+(R[3]-R[0])*b+(R[6]-R[0])*c;
  Result[1]=R[1]+(R[4]-R[1])*b+(R[7]-R[1])*c;
  Result[2]=R[2]+(R[5]-R[2])*b+(R[8]-R[2])*c;
  return;
}

template <typename ValueType>
void Interpolation_3D_Get_Volume(
ValueType * Result, 
const int N, 
const ValueType * R, 
const int I)  
{
  ValueType b=(ValueType)0;
  ValueType c=(ValueType)0;
  ValueType d=(ValueType)0;
  if(N==4) //0
  {
    b=(ValueType)0.25;//1
    c=(ValueType)0.25;//1
    d=(ValueType)0.25;//1
  }
  else if(N==5) //0,1,2,3
  {
    if(I==0)
    {
      b=(ValueType)0.2;//1
      c=(ValueType)0.2;//1
      d=(ValueType)0.2;//1
    }
    else if(I==1)
    {
      b=(ValueType)0.4;//2
      c=(ValueType)0.2;//1
      d=(ValueType)0.2;//1
    }
    else if(I==2)
    {
      b=(ValueType)0.2;//1
      c=(ValueType)0.4;//2
      d=(ValueType)0.2;//1
    }
    else
    {
      b=(ValueType)0.2;//1
      c=(ValueType)0.2;//1
      d=(ValueType)0.4;//2
    };
  };
  Result[0]=R[0]+(R[3]-R[0])*b+(R[6]-R[0])*c+(R[ 9]-R[0])*d;
  Result[1]=R[1]+(R[4]-R[1])*b+(R[7]-R[1])*c+(R[10]-R[1])*d;
  Result[2]=R[2]+(R[5]-R[2])*b+(R[8]-R[2])*c+(R[11]-R[2])*d;
  return;
}

//////////////////вспомогательные//////////////

template <typename ValueType>
ValueType interpolationVectorProduct2D(
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB)
{
  ValueType result;
  result=xA*yB-xB*yA;
  return(result);
}

template <typename ValueType>
ValueType interpolationCompositionalProduct3D(
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC)
{
  ValueType result;
  result=xA*yB*zC+xB*yC*zA+xC*yA*zB-xC*yB*zA-xB*yA*zC-xA*yC*zB;
  return(result);
}

template <typename ValueType>
void interpolationInitial_2D(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
ValueType * sA,
ValueType * sB,
ValueType * sC)
{
  ValueType S;
  S=interpolationVectorProduct2D<ValueType>(xB-xA, yB-yA, xC-xA, yC-yA);
  (*sA)=interpolationVectorProduct2D<ValueType>(xC-xB, yC-yB, x-xB, y-yB)/S;
  (*sB)=interpolationVectorProduct2D<ValueType>(xA-xC, yA-yC, x-xC, y-yC)/S;
  (*sC)=interpolationVectorProduct2D<ValueType>(xB-xA, yB-yA, x-xA, y-yA)/S;
  return;
}

template <typename ValueType>
void interpolationInitial_3D(
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
ValueType * sA,
ValueType * sB,
ValueType * sC,
ValueType * sD)
{
  ValueType S;
  S=interpolationCompositionalProduct3D<ValueType>(xC-xD, yC-yD, zC-zD, xA-xD, yA-yD, zA-zD, xB-xD, yB-yD, zB-zD);
  (*sA)=interpolationCompositionalProduct3D<ValueType>(x-xC, y-yC, z-zC, xD-xC, yD-yC, zD-zC, xB-xC, yB-yC, zB-zC)/S;
  (*sB)=interpolationCompositionalProduct3D<ValueType>(x-xD, y-yD, z-zD, xC-xD, yC-yD, zC-zD, xA-xD, yA-yD, zA-zD)/S;
  (*sC)=interpolationCompositionalProduct3D<ValueType>(x-xB, y-yB, z-zB, xD-xB, yD-yB, zD-zB, xA-xB, yA-yB, zA-zB)/S;
  (*sD)=interpolationCompositionalProduct3D<ValueType>(x-xC, y-yC, z-zC, xB-xC, yB-yC, zB-zC, xA-xC, yA-yC, zA-zC)/S;
  return;
}

//////////////////////////////////////////////N=1/////////////////////////////////////////

//2D

template <typename ValueType>
ValueType Interpolate_2D_Linear(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v100, const ValueType v010, const ValueType v001)
{
  ValueType result;
  ValueType sA, sB, sC;
  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);
  result=sA*v100+sB*v010+sC*v001;
  return(result);
}

//3D

template <typename ValueType>
ValueType Interpolate_3D_Linear(
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v1000, const ValueType v0100, const ValueType v0010, const ValueType v0001)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);
  result=sA*v1000+sB*v0100+sC*v0010+sD*v0001;
  return(result);
}

//////////////////////////////////////////////N=2//////////////////////////////////////////

//2D

template <typename ValueType>
ValueType Interpolate_2D_Power_2(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v200, const ValueType v020, const ValueType v002,
const ValueType v110, const ValueType v011, const ValueType v101)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType w200, w020, w002;
  ValueType w110, w011, w101;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);
  w200=sA*((ValueType)2.0*sA-(ValueType)1.0);
  w020=sB*((ValueType)2.0*sB-(ValueType)1.0);
  w002=sC*((ValueType)2.0*sC-(ValueType)1.0);

  w110=(ValueType)4.0*sA*sB;
  w011=(ValueType)4.0*sB*sC;
  w101=(ValueType)4.0*sC*sA;

  result=v200*w200+v020*w020+v002*w002+
         v110*w110+v011*w011+v101*w101;  
  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Linear_2(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v200, const ValueType v020, const ValueType v002,
const ValueType v110, const ValueType v011, const ValueType v101)
{
  ValueType result;
  ValueType sA, sB, sC;
  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);
  sA=(ValueType)2.0*sA;
  sB=(ValueType)2.0*sB;
  sC=(ValueType)2.0*sC;
  if (sA>(ValueType)1) //1
  {
    result=(sA-(ValueType)1.0)*v200+sB*v110+sC*v101;
  }
  else if (sB>(ValueType)1) //2
  {
    result=sA*v110+(sB-(ValueType)1.0)*v020+sC*v011;
  }
  else if (sC>(ValueType)1) //3
  {
    result=sA*v101+sB*v011+(sC-(ValueType)1.0)*v002;
  }
  else //4
  {
    result=((ValueType)1.0-sA)*v011+((ValueType)1.0-sB)*v101+((ValueType)1.0-sC)*v110;
  };
  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Mono_2(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v200, const ValueType v020, const ValueType v002,
const ValueType v110, const ValueType v011, const ValueType v101)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType Vmin, Vmax;
  ValueType w200, w020, w002;
  ValueType w110, w011, w101;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);

  w200=sA*((ValueType)2.0*sA-(ValueType)1.0);
  w020=sB*((ValueType)2.0*sB-(ValueType)1.0);
  w002=sC*((ValueType)2.0*sC-(ValueType)1.0);

  w110=(ValueType)4.0*sA*sB;
  w011=(ValueType)4.0*sB*sC;
  w101=(ValueType)4.0*sC*sA;

  result=v200*w200+v020*w020+v002*w002+
         v110*w110+v011*w011+v101*w101;  

  sA=(ValueType)2.0*sA;
  sB=(ValueType)2.0*sB;
  sC=(ValueType)2.0*sC;
  if (sA>(ValueType)1) //1
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v200, v110, v101);  
  }
  else if (sB>(ValueType)1) //2
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v110, v020, v011);  
  }
  else if (sC>(ValueType)1) //3
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v101, v011, v002);  
  }
  else //4
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v011, v101, v110);  
  };
  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };
  return(result);
}

//3D

template <typename ValueType>
ValueType Interpolate_3D_Power_2(
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v2000, const ValueType v0200, const ValueType v0020, const ValueType v0002,
const ValueType v1100, const ValueType v0110, const ValueType v0011, 
const ValueType v1001, const ValueType v1010, const ValueType v0101)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType w2000, w0200, w0020, w0002;
  ValueType w1100, w0110, w0011, w1001, w1010, w0101;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);
  w2000=sA*((ValueType)2.0*sA-(ValueType)1.0);
  w0200=sB*((ValueType)2.0*sB-(ValueType)1.0);
  w0020=sC*((ValueType)2.0*sC-(ValueType)1.0);
  w0002=sD*((ValueType)2.0*sD-(ValueType)1.0);

  w1100=(ValueType)4.0*sA*sB;
  w0110=(ValueType)4.0*sB*sC;
  w0011=(ValueType)4.0*sC*sD;
  w1001=(ValueType)4.0*sD*sA;
  w1010=(ValueType)4.0*sA*sC;
  w0101=(ValueType)4.0*sB*sD;

  result=v2000*w2000+v0200*w0200+v0020*w0020+v0002*w0002+
         v1100*w1100+v0110*w0110+v0011*w0011+v1001*w1001+v1010*w1010+v0101*w0101;
  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Linear_2(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v2000, const ValueType v0200, const ValueType v0020, const ValueType v0002,
const ValueType v1100, const ValueType v0110, const ValueType v0011, 
const ValueType v1001, const ValueType v1010, const ValueType v0101)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);
  sA=(ValueType)2.0*sA;
  sB=(ValueType)2.0*sB;
  sC=(ValueType)2.0*sC;
  sD=(ValueType)2.0*sD;  

  if (sA>(ValueType)1) //1
  {
    result=(sA-(ValueType)1.0)*v2000+sB*v1100+sC*v1010+sD*v1001;
  }
  else if (sB>(ValueType)1) //2
  {
    result=sA*v1100+(sB-(ValueType)1.0)*v0200+sC*v0110+sD*v0101;
  }
  else if (sC>(ValueType)1) //3
  {
    result=sA*v1010+sB*v0110+(sC-(ValueType)1.0)*v0020+sD*v0011;
  }
  else if (sD>(ValueType)1) //4
  {
    result=sA*v1001+sB*v0101+sC*v0011+(sD-(ValueType)1.0)*v0002;
  }
  else //5
  {
    if(axis=='0') //5.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0110+sD*v1001+s1*v1100+s2*v1010; //5.1.1
        }
        else
        {
          result=sC*v0110+((ValueType)1.0-sB)*v1001+s1*v1100-s2*v0101; //5.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0110+((ValueType)1.0-sC)*v1001-s1*v0011+s2*v1010; //5.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0110+sA*v1001-s1*v0011-s2*v0101; //5.1.2
        };
      };
    }
    else if(axis=='1') //5.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1010+sD*v0101+s1*v1100+s2*v0110; //5.2.1
        }
        else
        {
          result=sC*v1010+((ValueType)1.0-sA)*v0101+s1*v1100-s2*v1001; //5.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1010+((ValueType)1.0-sC)*v0101-s1*v0011+s2*v0110; //5.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1010+sB*v0101-s1*v0011-s2*v1001; //5.2.4
        };
      };
    }
    else //5.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1100+sD*v0011+s1*v0110+s2*v1010; //5.3.1
        }
        else
        {
          result=sA*v1100+((ValueType)1.0-sB)*v0011+s1*v0110-s2*v0101; //5.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1100+((ValueType)1.0-sA)*v0011-s1*v1001+s2*v1010; //5.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1100+sC*v0011-s1*v1001-s2*v0101; //5.3.3
        };
      };
    };        
  };
  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Mono_2(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v2000, const ValueType v0200, const ValueType v0020, const ValueType v0002,
const ValueType v1100, const ValueType v0110, const ValueType v0011, 
const ValueType v1001, const ValueType v1010, const ValueType v0101)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;
  ValueType Vmin, Vmax;
  ValueType w2000, w0200, w0020, w0002;
  ValueType w1100, w0110, w0011, w1001, w1010, w0101;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  w2000=sA*((ValueType)2.0*sA-(ValueType)1.0);
  w0200=sB*((ValueType)2.0*sB-(ValueType)1.0);
  w0020=sC*((ValueType)2.0*sC-(ValueType)1.0);
  w0002=sD*((ValueType)2.0*sD-(ValueType)1.0);

  w1100=(ValueType)4.0*sA*sB;
  w0110=(ValueType)4.0*sB*sC;
  w0011=(ValueType)4.0*sC*sD;
  w1001=(ValueType)4.0*sD*sA;
  w1010=(ValueType)4.0*sA*sC;
  w0101=(ValueType)4.0*sB*sD;

  result=v2000*w2000+v0200*w0200+v0020*w0020+v0002*w0002+
         v1100*w1100+v0110*w0110+v0011*w0011+v1001*w1001+v1010*w1010+v0101*w0101;

  sA=(ValueType)2.0*sA;
  sB=(ValueType)2.0*sB;
  sC=(ValueType)2.0*sC;
  sD=(ValueType)2.0*sD;  

  if (sA>(ValueType)1) //1
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2000, v1100, v1010, v1001);
  }
  else if (sB>(ValueType)1) //2
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1100, v0200, v0110, v0101);
  }
  else if (sC>(ValueType)1) //3
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1010, v0110, v0020, v0011);
  }
  else if (sD>(ValueType)1) //4
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1001, v0101, v0011, v0002);
  }
  else //5
  {
    if(axis=='0') //5.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0110, v1001, v1100, v1010); //5.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0110, v1001, v1100, v0101); //5.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0110, v1001, v0011, v1010); //5.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0110, v1001, v0011, v0101); //5.1.2
        };
      };
    }
    else if(axis=='1') //5.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1010, v0101, v1100, v0110); //5.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1010, v0101, v1100, v1001); //5.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1010, v0101, v0011, v0110); //5.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1010, v0101, v0011, v1001); //5.2.4
        };
      };
    }
    else //5.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1100, v0011, v0110, v1010); //5.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1100, v0011, v0110, v0101); //5.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1100, v0011, v1001, v1010); //5.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1100, v0011, v1001, v0101); //5.3.3
        };
      };
    };        
  };

  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };

  return(result);
}

///////////////////////////////////////////////N=3//////////////////////////////////////////

//2D

template <typename ValueType>
ValueType Interpolate_2D_Power_3(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v300, const ValueType v030, const ValueType v003,
const ValueType v210, const ValueType v021, const ValueType v102, 
const ValueType v120, const ValueType v012, const ValueType v201,
const ValueType v111)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType w300, w030, w003;
  ValueType w210, w021, w102; 
  ValueType w120, w012, w201;
  ValueType w111;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);

  w300=sA*((ValueType)3.0*sA-(ValueType)1.0)*((ValueType)3.0*sA-(ValueType)2.0)/(ValueType)2.0;
  w030=sB*((ValueType)3.0*sB-(ValueType)1.0)*((ValueType)3.0*sB-(ValueType)2.0)/(ValueType)2.0;
  w003=sC*((ValueType)3.0*sC-(ValueType)1.0)*((ValueType)3.0*sC-(ValueType)2.0)/(ValueType)2.0;

  w210=sA*sB*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w021=sB*sC*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w102=sC*sA*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w120=sB*sA*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w012=sC*sB*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w201=sA*sC*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);

  w111=(ValueType)27.0*sA*sB*sC;

  result=v300*w300+v030*w030+v003*w003+
         v210*w210+v021*w021+v102*w102+
         v120*w120+v012*w012+v201*w201+
         v111*w111;
  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Linear_3(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v300, const ValueType v030, const ValueType v003,
const ValueType v210, const ValueType v021, const ValueType v102, 
const ValueType v120, const ValueType v012, const ValueType v201,
const ValueType v111)
{
  ValueType result;
  ValueType sA, sB, sC;
  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);
  sA=(ValueType)3.0*sA;
  sB=(ValueType)3.0*sB;
  sC=(ValueType)3.0*sC;

  if(sA>(ValueType)2) //1
  {
    result=(sA-(ValueType)2.0)*v300+sB*v210+sC*v201;
  }
  else if(sB>(ValueType)2) //2
  {
    result=sA*v120+(sB-(ValueType)2.0)*v030+sC*v021;
  }
  else if(sC>(ValueType)2) //3
  {
    result=sA*v102+sB*v012+(sC-(ValueType)2.0)*v003;
  }
  else if ((sB<(ValueType)1)&&(sC<(ValueType)1)) //4
  {
    result=((ValueType)2.0-sA)*v111+((ValueType)1.0-sB)*v201+((ValueType)1.0-sC)*v210;
  }
  else if ((sC<(ValueType)1)&&(sA<(ValueType)1)) //5
  {
    result=((ValueType)1.0-sA)*v021+((ValueType)2.0-sB)*v111+((ValueType)1.0-sC)*v120;
  }
  else if ((sA<(ValueType)1)&&(sB<(ValueType)1)) //6
  {
    result=((ValueType)1.0-sA)*v012+((ValueType)1.0-sB)*v102+((ValueType)2.0-sC)*v111;
  }
  else if (sA<(ValueType)1) //7
  {
    result=sA*v111+(sB-(ValueType)1.0)*v021+(sC-(ValueType)1.0)*v012;
  }
  else if (sB<(ValueType)1) //8
  {
    result=(sA-(ValueType)1.0)*v201+sB*v111+(sC-(ValueType)1.0)*v102;
  }
  else //9
  {
    result=(sA-(ValueType)1.0)*v210+(sB-(ValueType)1.0)*v120+sC*v111;
  };
  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Mono_3(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v300, const ValueType v030, const ValueType v003,
const ValueType v210, const ValueType v021, const ValueType v102, 
const ValueType v120, const ValueType v012, const ValueType v201,
const ValueType v111)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType Vmin, Vmax;
  ValueType w300, w030, w003;
  ValueType w210, w021, w102; 
  ValueType w120, w012, w201;
  ValueType w111;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);
  w300=sA*((ValueType)3.0*sA-(ValueType)1.0)*((ValueType)3.0*sA-(ValueType)2.0)/(ValueType)2.0;
  w030=sB*((ValueType)3.0*sB-(ValueType)1.0)*((ValueType)3.0*sB-(ValueType)2.0)/(ValueType)2.0;
  w003=sC*((ValueType)3.0*sC-(ValueType)1.0)*((ValueType)3.0*sC-(ValueType)2.0)/(ValueType)2.0;

  w210=sA*sB*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w021=sB*sC*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w102=sC*sA*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w120=sB*sA*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w012=sC*sB*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w201=sA*sC*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);

  w111=(ValueType)27.0*sA*sB*sC;

  result=v300*w300+v030*w030+v003*w003+
         v210*w210+v021*w021+v102*w102+
         v120*w120+v012*w012+v201*w201+
         v111*w111;

  sA=(ValueType)3.0*sA;
  sB=(ValueType)3.0*sB;
  sC=(ValueType)3.0*sC;

  if(sA>(ValueType)2) //1
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v300, v210, v201);  
  }
  else if(sB>(ValueType)2) //2
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v120, v030, v021);  
  }
  else if(sC>(ValueType)2) //3
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v102, v012, v003);  
  }
  else if ((sB<(ValueType)1)&&(sC<(ValueType)1)) //4
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v111, v201, v210);  
  }
  else if ((sC<(ValueType)1)&&(sA<(ValueType)1)) //5
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v021, v111, v120);  
  }
  else if ((sA<(ValueType)1)&&(sB<(ValueType)1)) //6
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v012, v102, v111);  
  }
  else if (sA<(ValueType)1) //7
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v111, v021, v012);  
  }
  else if (sB<(ValueType)1) //8
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v201, v111, v102);  
  }
  else //9
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v210, v120, v111);  
  };
  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };
  return(result);
  return(result);
}

//3D

template <typename ValueType>
ValueType Interpolate_3D_Power_3(
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v3000, const ValueType v0300, const ValueType v0030, const ValueType v0003,
const ValueType v2100, const ValueType v0210, const ValueType v0021, const ValueType v1002, 
const ValueType v1200, const ValueType v0120, const ValueType v0012, const ValueType v2001, 
const ValueType v2010, const ValueType v0201, const ValueType v1020, const ValueType v0102,
const ValueType v1110, const ValueType v0111, const ValueType v1011, const ValueType v1101)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType w3000, w0300, w0030, w0003;
  ValueType w2100, w0210, w0021, w1002; 
  ValueType w1200, w0120, w0012, w2001; 
  ValueType w2010, w0201, w1020, w0102;
  ValueType w1110, w0111, w1011, w1101;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  w3000=sA*((ValueType)3.0*sA-(ValueType)1.0)*((ValueType)3.0*sA-(ValueType)2.0)/(ValueType)2.0;
  w0300=sB*((ValueType)3.0*sB-(ValueType)1.0)*((ValueType)3.0*sB-(ValueType)2.0)/(ValueType)2.0;
  w0030=sC*((ValueType)3.0*sC-(ValueType)1.0)*((ValueType)3.0*sC-(ValueType)2.0)/(ValueType)2.0;
  w0003=sD*((ValueType)3.0*sD-(ValueType)1.0)*((ValueType)3.0*sD-(ValueType)2.0)/(ValueType)2.0;

  w2100=sA*sB*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0210=sB*sC*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0021=sC*sD*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w1002=sD*sA*((ValueType)3.0*sD-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w1200=sB*sA*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0120=sC*sB*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0012=sD*sC*((ValueType)3.0*sD-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w2001=sA*sD*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w2010=sA*sC*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0201=sB*sD*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w1020=sC*sA*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0102=sD*sB*((ValueType)3.0*sD-(ValueType)1.0)*(ValueType)(9.0/2.0);

  w1110=(ValueType)27.0*sA*sB*sC;
  w0111=(ValueType)27.0*sB*sC*sD;
  w1011=(ValueType)27.0*sC*sD*sA;
  w1101=(ValueType)27.0*sD*sA*sB;

  result=v3000*w3000+v0300*w0300+v0030*w0030+v0003*w0003+
         v2100*w2100+v0210*w0210+v0021*w0021+v1002*w1002+
         v1200*w1200+v0120*w0120+v0012*w0012+v2001*w2001+
         v2010*w2010+v0201*w0201+v1020*w1020+v0102*w0102+
         v1110*w1110+v0111*w0111+v1011*w1011+v1101*w1101;
  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Linear_3(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v3000, const ValueType v0300, const ValueType v0030, const ValueType v0003,
const ValueType v2100, const ValueType v0210, const ValueType v0021, const ValueType v1002, 
const ValueType v1200, const ValueType v0120, const ValueType v0012, const ValueType v2001, 
const ValueType v2010, const ValueType v0201, const ValueType v1020, const ValueType v0102,
const ValueType v1110, const ValueType v0111, const ValueType v1011, const ValueType v1101)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);
  sA=(ValueType)3.0*sA;
  sB=(ValueType)3.0*sB;
  sC=(ValueType)3.0*sC;
  sD=(ValueType)3.0*sD;  

  if (sA>(ValueType)2) //1
  {
    result=(sA-(ValueType)2.0)*v3000+sB*v2100+sC*v2010+sD*v2001;
  }
  else if(sB>(ValueType)2) //2
  {
    result=sA*v1200+(sB-(ValueType)2.0)*v0300+sC*v0210+sD*v0201;
  }
  else if(sC>(ValueType)2) //3
  {
    result=sA*v1020+sB*v0120+(sC-(ValueType)2.0)*v0030+sD*v0021;
  }
  else if(sD>(ValueType)2) //4
  {
    result=sA*v1002+sB*v0102+sC*v0012+(sD-(ValueType)2.0)*v0003;
  }
  else if((sA>(ValueType)1)&&(sB>(ValueType)1)) //5
  {
    result=(sA-(ValueType)1.0)*v2100+(sB-(ValueType)1.0)*v1200+sC*v1110+sD*v1101;
  }
  else if((sB>(ValueType)1)&&(sC>(ValueType)1)) //6
  {
    result=sA*v1110+(sB-(ValueType)1.0)*v0210+(sC-(ValueType)1.0)*v0120+sD*v0111;
  }
  else if((sC>(ValueType)1)&&(sD>(ValueType)1)) //7
  {
    result=sA*v1011+sB*v0111+(sC-(ValueType)1.0)*v0021+(sD-(ValueType)1.0)*v0012;
  }
  else if((sA>(ValueType)1)&&(sD>(ValueType)1)) //8
  {
    result=(sA-(ValueType)1.0)*v2001+sB*v1101+sC*v1011+(sD-(ValueType)1.0)*v1002;
  }
  else if((sA>(ValueType)1)&&(sC>(ValueType)1)) //9
  {
    result=(sA-(ValueType)1.0)*v2010+sB*v1110+(sC-(ValueType)1.0)*v1020+sD*v1011;
  }
  else if((sB>(ValueType)1)&&(sD>(ValueType)1)) //10
  {
    result=sA*v1101+(sB-(ValueType)1.0)*v0201+sC*v0111+(sD-(ValueType)1.0)*v0102;
  }
  else if(sA>(ValueType)1) //11
  {
    if(axis=='0') //11.1
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1110+sD*v2001+s1*v2100+s2*v2010; //11.1.1
        }
        else
        {
          result=sC*v1110+((ValueType)1.0-sB)*v2001+s1*v2100-s2*v1101; //11.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1110+((ValueType)1.0-sC)*v2001-s1*v1011+s2*v2010; //11.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1110+(sA-(ValueType)1.0)*v2001-s1*v1011-s2*v1101; //11.1.2
        };
      };
    }
    else if(axis=='1') //11.2
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v2010+sD*v1101+s1*v2100+s2*v1110; //11.2.1
        }
        else
        {
          result=sC*v2010+((ValueType)2.0-sA)*v1101+s1*v2100-s2*v2001; //11.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2010+((ValueType)1.0-sC)*v1101-s1*v1011+s2*v1110; //11.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v2010+sB*v1101-s1*v1011-s2*v2001; //11.2.4
        };
      };
    }
    else //11.3
    {
      s1=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v2100+sD*v1011+s1*v1110+s2*v2010; //11.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2100+((ValueType)1.0-sB)*v1011+s1*v1110-s2*v1101; //11.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2100+((ValueType)2.0-sA)*v1011-s1*v2001+s2*v2010; //11.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2100+sC*v1011-s1*v2001-s2*v1101; //11.3.3
        };
      };
    };     
  }
  else if(sB>(ValueType)1) //12
  {
    if(axis=='0') //12.1
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0210+sD*v1101+s1*v1200+s2*v1110; //12.1.1
        }
        else
        {
          result=sC*v0210+((ValueType)2.0-sB)*v1101+s1*v1200-s2*v0201; //12.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v0210+((ValueType)1.0-sC)*v1101-s1*v0111+s2*v1110; //12.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0210+sA*v1101-s1*v0111-s2*v0201; //12.1.2
        };
      };
    }
    else if(axis=='1') //12.2
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v1110+sD*v0201+s1*v1200+s2*v0210; //12.2.1
        }
        else
        {
          result=sC*v1110+((ValueType)1.0-sA)*v0201+s1*v1200-s2*v1101; //12.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1110+((ValueType)1.0-sC)*v0201-s1*v0111+s2*v0210; //12.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1110+(sB-(ValueType)1.0)*v0201-s1*v0111-s2*v1101; //12.2.4
        };
      };
    }
    else //12.3
    {
      s1=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1200+sD*v0111+s1*v0210+s2*v1110; //12.3.1
        }
        else
        {
          result=sA*v1200+((ValueType)2.0-sB)*v0111+s1*v0210-s2*v0201; //12.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1200+((ValueType)1.0-sA)*v0111-s1*v1101+s2*v1110; //12.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1200+sC*v0111-s1*v1101-s2*v0201; //12.3.3
        };
      };
    };  
  }
  else if(sC>(ValueType)1) //13
  {
    if(axis=='0') //13.1
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0120+sD*v1011+s1*v1110+s2*v1020; //13.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v0120+((ValueType)1.0-sB)*v1011+s1*v1110-s2*v0111; //13.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0120+((ValueType)2.0-sC)*v1011-s1*v0021+s2*v1020; //13.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0120+sA*v1011-s1*v0021-s2*v0111; //13.1.2
        };
      };
    }
    else if(axis=='1') //13.2
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1020+sD*v0111+s1*v1110+s2*v0120; //13.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1020+((ValueType)1.0-sA)*v0111+s1*v1110-s2*v1011; //13.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1020+((ValueType)2.0-sC)*v0111-s1*v0021+s2*v0120; //13.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1020+sB*v0111-s1*v0021-s2*v1011; //13.2.4
        };
      };
    }
    else //13.3
    {
      s1=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v1110+sD*v0021+s1*v0120+s2*v1020; //13.3.1
        }
        else
        {
          result=sA*v1110+((ValueType)1.0-sB)*v0021+s1*v0120-s2*v0111; //13.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1110+((ValueType)1.0-sA)*v0021-s1*v1011+s2*v1020; //13.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1110+(sC-(ValueType)1.0)*v0021-s1*v1011-s2*v0111; //13.3.3
        };
      };
    };  
  }      
  else if(sD>(ValueType)1)//14
  {
    if(axis=='0') //14.1
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0111+(sD-(ValueType)1.0)*v1002+s1*v1101+s2*v1011; //14.1.1
        }
        else
        {
          result=sC*v0111+((ValueType)1.0-sB)*v1002+s1*v1101-s2*v0102; //14.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0111+((ValueType)1.0-sC)*v1002-s1*v0012+s2*v1011; //14.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v0111+sA*v1002-s1*v0012-s2*v0102; //14.1.2
        };
      };
    }
    else if(axis=='1') //14.2
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1011+(sD-(ValueType)1.0)*v0102+s1*v1101+s2*v0111; //14.2.1
        }
        else
        {
          result=sC*v1011+((ValueType)1.0-sA)*v0102+s1*v1101-s2*v1002; //14.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1011+((ValueType)1.0-sC)*v0102-s1*v0012+s2*v0111; //14.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v1011+sB*v0102-s1*v0012-s2*v1002; //14.2.4
        };
      };
    }
    else //14.3
    {
      s1=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1101+(sD-(ValueType)1.0)*v0012+s1*v0111+s2*v1011; //14.3.1
        }
        else
        {
          result=sA*v1101+((ValueType)1.0-sB)*v0012+s1*v0111-s2*v0102; //14.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1101+((ValueType)1.0-sA)*v0012-s1*v1002+s2*v1011; //14.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1101+sC*v0012-s1*v1002-s2*v0102; //14.3.3
        };
      };
    };
  }  
  else //15
  {
    result=((ValueType)1-sA)*v0111+((ValueType)1-sB)*v1011+((ValueType)1-sC)*v1101+((ValueType)1-sD)*v1110;
  };
  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Mono_3(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v3000, const ValueType v0300, const ValueType v0030, const ValueType v0003,
const ValueType v2100, const ValueType v0210, const ValueType v0021, const ValueType v1002, 
const ValueType v1200, const ValueType v0120, const ValueType v0012, const ValueType v2001, 
const ValueType v2010, const ValueType v0201, const ValueType v1020, const ValueType v0102,
const ValueType v1110, const ValueType v0111, const ValueType v1011, const ValueType v1101)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;
  ValueType Vmin, Vmax;
  ValueType w3000, w0300, w0030, w0003;
  ValueType w2100, w0210, w0021, w1002; 
  ValueType w1200, w0120, w0012, w2001; 
  ValueType w2010, w0201, w1020, w0102;
  ValueType w1110, w0111, w1011, w1101;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  w3000=sA*((ValueType)3.0*sA-(ValueType)1.0)*((ValueType)3.0*sA-(ValueType)2.0)/(ValueType)2.0;
  w0300=sB*((ValueType)3.0*sB-(ValueType)1.0)*((ValueType)3.0*sB-(ValueType)2.0)/(ValueType)2.0;
  w0030=sC*((ValueType)3.0*sC-(ValueType)1.0)*((ValueType)3.0*sC-(ValueType)2.0)/(ValueType)2.0;
  w0003=sD*((ValueType)3.0*sD-(ValueType)1.0)*((ValueType)3.0*sD-(ValueType)2.0)/(ValueType)2.0;

  w2100=sA*sB*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0210=sB*sC*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0021=sC*sD*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w1002=sD*sA*((ValueType)3.0*sD-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w1200=sB*sA*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0120=sC*sB*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0012=sD*sC*((ValueType)3.0*sD-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w2001=sA*sD*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w2010=sA*sC*((ValueType)3.0*sA-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0201=sB*sD*((ValueType)3.0*sB-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w1020=sC*sA*((ValueType)3.0*sC-(ValueType)1.0)*(ValueType)(9.0/2.0);
  w0102=sD*sB*((ValueType)3.0*sD-(ValueType)1.0)*(ValueType)(9.0/2.0);

  w1110=(ValueType)27.0*sA*sB*sC;
  w0111=(ValueType)27.0*sB*sC*sD;
  w1011=(ValueType)27.0*sC*sD*sA;
  w1101=(ValueType)27.0*sD*sA*sB;

  result=v3000*w3000+v0300*w0300+v0030*w0030+v0003*w0003+
         v2100*w2100+v0210*w0210+v0021*w0021+v1002*w1002+
         v1200*w1200+v0120*w0120+v0012*w0012+v2001*w2001+
         v2010*w2010+v0201*w0201+v1020*w1020+v0102*w0102+
         v1110*w1110+v0111*w0111+v1011*w1011+v1101*w1101;

  sA=(ValueType)3.0*sA;
  sB=(ValueType)3.0*sB;
  sC=(ValueType)3.0*sC;
  sD=(ValueType)3.0*sD;  

  if (sA>(ValueType)2) //1
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3000, v2100, v2010, v2001);
  }
  else if(sB>(ValueType)2) //2
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1200, v0300, v0210, v0201);
  }
  else if(sC>(ValueType)2) //3
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1020, v0120, v0030, v0021);
  }
  else if(sD>(ValueType)2) //4
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1002, v0102, v0012, v0003);
  }
  else if((sA>(ValueType)1)&&(sB>(ValueType)1)) //5
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2100, v1200, v1110, v1101);
  }
  else if((sB>(ValueType)1)&&(sC>(ValueType)1)) //6
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0210, v0120, v0111);
  }
  else if((sC>(ValueType)1)&&(sD>(ValueType)1)) //7
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1011, v0111, v0021, v0012);
  }
  else if((sA>(ValueType)1)&&(sD>(ValueType)1)) //8
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2001, v1101, v1011, v1002);
  }
  else if((sA>(ValueType)1)&&(sC>(ValueType)1)) //9
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2010, v1110, v1020, v1011);
  }
  else if((sB>(ValueType)1)&&(sD>(ValueType)1)) //10
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1101, v0201, v0111, v0102);
  }
  else if(sA>(ValueType)1) //11
  {
    if(axis=='0') //11.1
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v2001, v2100, v2010); //11.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v2001, v2100, v1101); //11.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v2001, v1011, v2010); //11.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v2001, v1011, v1101); //11.1.2
        };
      };
    }
    else if(axis=='1') //11.2
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2010, v1101, v2100, v1110); //11.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2010, v1101, v2100, v2001); //11.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2010, v1101, v1011, v1110); //11.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2010, v1101, v1011, v2001); //11.2.4
        };
      };
    }
    else //11.3
    {
      s1=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2100, v1011, v1110, v2010); //11.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2100, v1011, v1110, v1101); //11.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2100, v1011, v2001, v2010); //11.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2100, v1011, v2001, v1101); //11.3.3
        };
      };
    };     
  }
  else if(sB>(ValueType)1) //12
  {
    if(axis=='0') //12.1
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0210, v1101, v1200, v1110); //12.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0210, v1101, v1200, v0201); //12.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0210, v1101, v0111, v1110); //12.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0210, v1101, v0111, v0201); //12.1.2
        };
      };
    }
    else if(axis=='1') //12.2
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0201, v1200, v0210); //12.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0201, v1200, v1101); //12.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0201, v0111, v0210); //12.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0201, v0111, v1101); //12.2.4
        };
      };
    }
    else //12.3
    {
      s1=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1200, v0111, v0210, v1110); //12.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1200, v0111, v0210, v0201); //12.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1200, v0111, v1101, v1110); //12.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1200, v0111, v1101, v0201); //12.3.3
        };
      };
    };  
  }
  else if(sC>(ValueType)1) //13
  {
    if(axis=='0') //13.1
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0120, v1011, v1110, v1020); //13.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0120, v1011, v1110, v0111); //13.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0120, v1011, v0021, v1020); //13.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0120, v1011, v0021, v0111); //13.1.2
        };
      };
    }
    else if(axis=='1') //13.2
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1020, v0111, v1110, v0120); //13.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1020, v0111, v1110, v1011); //13.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1020, v0111, v0021, v0120); //13.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1020, v0111, v0021, v1011); //13.2.4
        };
      };
    }
    else //13.3
    {
      s1=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0021, v0120, v1020); //13.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0021, v0120, v0111); //13.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0021, v1011, v1020); //13.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1110, v0021, v1011, v0111); //13.3.3
        };
      };
    };  
  }      
  else if(sD>(ValueType)1)//14
  {
    if(axis=='0') //14.1
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0111, v1002, v1101, v1011); //14.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0111, v1002, v1101, v0102); //14.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0111, v1002, v0012, v1011); //14.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0111, v1002, v0012, v0102); //14.1.2
        };
      };
    }
    else if(axis=='1') //14.2
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1011, v0102, v1101, v0111); //14.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1011, v0102, v1101, v1002); //14.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1011, v0102, v0012, v0111); //14.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1011, v0102, v0012, v1002); //14.2.4
        };
      };
    }
    else //14.3
    {
      s1=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1101, v0012, v0111, v1011); //14.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1101, v0012, v0111, v0102); //14.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1101, v0012, v1002, v1011); //14.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1101, v0012, v1002, v0102); //14.3.3
        };
      };
    };
  }  
  else //15
  {
    interpolationMinMax_3D(&Vmin, &Vmax, v0111, v1011, v1101, v1110);
  };

  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };

  return(result);
}

///////////////////////////////////////////////N=4//////////////////////////////////////////

//2D

template <typename ValueType>
ValueType Interpolate_2D_Power_4(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v400, const ValueType v040, const ValueType v004,
const ValueType v310, const ValueType v031, const ValueType v103, 
const ValueType v130, const ValueType v013, const ValueType v301,
const ValueType v220, const ValueType v022, const ValueType v202,
const ValueType v211, const ValueType v121, const ValueType v112)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType w400, w040, w004;
  ValueType w310, w031, w103; 
  ValueType w130, w013, w301;
  ValueType w220, w022, w202;
  ValueType w211, w121, w112;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);

  w400=sA*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)3.0)/(ValueType)3.0;
  w040=sB*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)3.0)/(ValueType)3.0;
  w004=sC*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)3.0)/(ValueType)3.0;

  w310=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w031=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w103=sC*sA*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w130=sB*sA*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w013=sC*sB*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w301=sA*sC*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);

  w220=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)4.0;
  w022=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)4.0;
  w202=sC*sA*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)4.0;

  w211=sA*sB*sC*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w121=sB*sC*sA*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;
  w112=sC*sA*sB*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;

  result=v400*w400+v040*w040+v004*w004+
         v310*w310+v031*w031+v103*w103+
         v130*w130+v013*w013+v301*w301+
         v220*w220+v022*w022+v202*w202+
         v211*w211+v121*w121+v112*w112;
  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Linear_4(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v400, const ValueType v040, const ValueType v004,
const ValueType v310, const ValueType v031, const ValueType v103, 
const ValueType v130, const ValueType v013, const ValueType v301,
const ValueType v220, const ValueType v022, const ValueType v202,
const ValueType v211, const ValueType v121, const ValueType v112)
{
  ValueType result;
  ValueType sA, sB, sC;
  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);
  sA=(ValueType)4.0*sA;
  sB=(ValueType)4.0*sB;
  sC=(ValueType)4.0*sC;

  if (sA>(ValueType)3) //1
  {
    result=(sA-(ValueType)3.0)*v400+sB*v310+sC*v301;
  }
  else if(sB>(ValueType)3) //2
  {
    result=sA*v130+(sB-(ValueType)3.0)*v040+sC*v031;
  }
  else if(sC>(ValueType)3) //3
  {
    result=sA*v103+sB*v013+(sC-(ValueType)3.0)*v004;
  }
  else if((sB<(ValueType)1)&&(sC<(ValueType)1)) //4
  {
    result=((ValueType)3.0-sA)*v211+((ValueType)1.0-sB)*v301+((ValueType)1.0-sC)*v310;
  }
  else if((sC<(ValueType)1)&&(sA<(ValueType)1)) //5
  {
    result=((ValueType)1.0-sA)*v031+((ValueType)3.0-sB)*v121+((ValueType)1.0-sC)*v130;
  }
  else if((sA<(ValueType)1)&&(sB<(ValueType)1)) //6
  {
    result=((ValueType)1.0-sA)*v013+((ValueType)1.0-sB)*v103+((ValueType)3.0-sC)*v112;
  }
  else if((sA<(ValueType)1)&&(sC>(ValueType)2)) //7
  {
    result=sA*v112+(sB-(ValueType)1.0)*v022+(sC-(ValueType)2.0)*v013;
  }
  else if((sB<(ValueType)1)&&(sC>(ValueType)2)) //8
  {
    result=(sA-(ValueType)1.0)*v202+sB*v112+(sC-(ValueType)2.0)*v103;
  }
  else if((sB<(ValueType)1)&&(sA>(ValueType)2)) //9
  {
    result=(sA-(ValueType)2.0)*v301+sB*v211+(sC-(ValueType)1.0)*v202;
  }
  else if((sC<(ValueType)1)&&(sA>(ValueType)2)) //10
  {
    result=(sA-(ValueType)2.0)*v310+(sB-(ValueType)1.0)*v220+sC*v211;
  }
  else if((sC<(ValueType)1)&&(sB>(ValueType)2)) //11
  {
    result=(sA-(ValueType)1.0)*v220+(sB-(ValueType)2.0)*v130+sC*v121;
  }
  else if((sA<(ValueType)1)&&(sB>(ValueType)2)) //12
  {
    result=sA*v121+(sB-(ValueType)2.0)*v031+(sC-(ValueType)1.0)*v022;
  }
  else if(sA<(ValueType)1) //13
  {
    result=((ValueType)1.0-sA)*v022+((ValueType)2.0-sB)*v112+((ValueType)2.0-sC)*v121;
  }
  else if(sB<(ValueType)1) //14
  {
    result=((ValueType)2.0-sA)*v112+((ValueType)1.0-sB)*v202+((ValueType)2.0-sC)*v211;
  }
  else if(sC<(ValueType)1) //15
  {
    result=((ValueType)2.0-sA)*v121+((ValueType)2.0-sB)*v211+((ValueType)1.0-sC)*v220;
  }
  else //16
  {
    result=(sA-(ValueType)1.0)*v211+(sB-(ValueType)1.0)*v121+(sC-(ValueType)1.0)*v112;
  };
  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Mono_4(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v400, const ValueType v040, const ValueType v004,
const ValueType v310, const ValueType v031, const ValueType v103, 
const ValueType v130, const ValueType v013, const ValueType v301,
const ValueType v220, const ValueType v022, const ValueType v202,
const ValueType v211, const ValueType v121, const ValueType v112)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType Vmin, Vmax;
  ValueType w400, w040, w004;
  ValueType w310, w031, w103; 
  ValueType w130, w013, w301;
  ValueType w220, w022, w202;
  ValueType w211, w121, w112;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);

  w400=sA*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)3.0)/(ValueType)3.0;
  w040=sB*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)3.0)/(ValueType)3.0;
  w004=sC*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)3.0)/(ValueType)3.0;

  w310=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w031=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w103=sC*sA*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w130=sB*sA*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w013=sC*sB*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w301=sA*sC*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);

  w220=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)4.0;
  w022=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)4.0;
  w202=sC*sA*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)4.0;

  w211=sA*sB*sC*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w121=sB*sC*sA*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;
  w112=sC*sA*sB*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;

  result=v400*w400+v040*w040+v004*w004+
         v310*w310+v031*w031+v103*w103+
         v130*w130+v013*w013+v301*w301+
         v220*w220+v022*w022+v202*w202+
         v211*w211+v121*w121+v112*w112;

  sA=(ValueType)4.0*sA;
  sB=(ValueType)4.0*sB;
  sC=(ValueType)4.0*sC;

  if (sA>(ValueType)3) //1
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v400, v310, v301);  
  }
  else if(sB>(ValueType)3) //2
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v130, v040, v031);  
  }
  else if(sC>(ValueType)3) //3
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v103, v013, v004);  
  }
  else if((sB<(ValueType)1)&&(sC<(ValueType)1)) //4
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v211, v301, v310);  
  }
  else if((sC<(ValueType)1)&&(sA<(ValueType)1)) //5
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v031, v121, v130);  
  }
  else if((sA<(ValueType)1)&&(sB<(ValueType)1)) //6
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v013, v103, v112);  
  }
  else if((sA<(ValueType)1)&&(sC>(ValueType)2)) //7
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v112, v022, v013);  
  }
  else if((sB<(ValueType)1)&&(sC>(ValueType)2)) //8
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v202, v112, v103);  
  }
  else if((sB<(ValueType)1)&&(sA>(ValueType)2)) //9
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v301, v211, v202);  
  }
  else if((sC<(ValueType)1)&&(sA>(ValueType)2)) //10
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v310, v220, v211);  
  }
  else if((sC<(ValueType)1)&&(sB>(ValueType)2)) //11
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v220, v130, v121);  
  }
  else if((sA<(ValueType)1)&&(sB>(ValueType)2)) //12
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v121, v031, v022);  
  }
  else if(sA<(ValueType)1) //13
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v022, v112, v121);  
  }
  else if(sB<(ValueType)1) //14
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v112, v202, v211);  
  }
  else if(sC<(ValueType)1) //15
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v121, v211, v220);  
  }
  else //16
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v211, v121, v112);  
  };
  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };
  return(result);
}

//3D

template <typename ValueType>
ValueType Interpolate_3D_Power_4(
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v4000, const ValueType v0400, const ValueType v0040, const ValueType v0004,
const ValueType v3100, const ValueType v0310, const ValueType v0031, const ValueType v1003, 
const ValueType v1300, const ValueType v0130, const ValueType v0013, const ValueType v3001, 
const ValueType v3010, const ValueType v0301, const ValueType v1030, const ValueType v0103,
const ValueType v2200, const ValueType v0220, const ValueType v0022, 
const ValueType v2002, const ValueType v2020, const ValueType v0202,
const ValueType v2110, const ValueType v0211, const ValueType v1021, const ValueType v1102, 
const ValueType v1210, const ValueType v0121, const ValueType v1012, const ValueType v2101, 
const ValueType v1120, const ValueType v0112, const ValueType v2011, const ValueType v1201,
const ValueType v1111)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType w4000, w0400, w0040, w0004;
  ValueType w3100, w0310, w0031, w1003; 
  ValueType w1300, w0130, w0013, w3001; 
  ValueType w3010, w0301, w1030, w0103;
  ValueType w2200, w0220, w0022, w2002, w2020, w0202;
  ValueType w2110, w0211, w1021, w1102; 
  ValueType w1210, w0121, w1012, w2101; 
  ValueType w1120, w0112, w2011, w1201;
  ValueType w1111;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  w4000=sA*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)3.0)/(ValueType)3.0;
  w0400=sB*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)3.0)/(ValueType)3.0;
  w0040=sC*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)3.0)/(ValueType)3.0;
  w0004=sD*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*((ValueType)4.0*sD-(ValueType)3.0)/(ValueType)3.0;

  w3100=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0310=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0031=sC*sD*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w1003=sD*sA*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w1300=sB*sA*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0130=sC*sB*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0013=sD*sC*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w3001=sA*sD*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w3010=sA*sC*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0301=sB*sD*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w1030=sC*sA*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0103=sD*sB*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*(ValueType)(16.0/3.0);

  w2200=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)4.0;
  w0220=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)4.0;
  w0022=sC*sD*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)4.0;
  w2002=sD*sA*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)4.0;
  w2020=sA*sC*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)4.0;
  w0202=sB*sD*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)4.0;

  w2110=sA*sB*sC*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w0211=sB*sC*sD*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;
  w1021=sC*sD*sA*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;
  w1102=sD*sA*sB*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)32.0;
  w1210=sB*sC*sA*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;
  w0121=sC*sD*sB*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;
  w1012=sD*sA*sC*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)32.0;
  w2101=sA*sB*sD*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w1120=sC*sA*sB*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;
  w0112=sD*sB*sC*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)32.0;
  w2011=sA*sC*sD*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w1201=sB*sD*sA*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;

  w1111=sA*sB*sC*sD*(ValueType)256.0;

  result=v4000*w4000+v0400*w0400+v0040*w0040+v0004*w0004+
         v3100*w3100+v0310*w0310+v0031*w0031+v1003*w1003+
         v1300*w1300+v0130*w0130+v0013*w0013+v3001*w3001+
         v3010*w3010+v0301*w0301+v1030*w1030+v0103*w0103+
         v2200*w2200+v0220*w0220+v0022*w0022+v2002*w2002+v2020*w2020+v0202*w0202+
         v2110*w2110+v0211*w0211+v1021*w1021+v1102*w1102+
         v1210*w1210+v0121*w0121+v1012*w1012+v2101*w2101+
         v1120*w1120+v0112*w0112+v2011*w2011+v1201*w1201+
         v1111*w1111;
  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Linear_4(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v4000, const ValueType v0400, const ValueType v0040, const ValueType v0004,
const ValueType v3100, const ValueType v0310, const ValueType v0031, const ValueType v1003, 
const ValueType v1300, const ValueType v0130, const ValueType v0013, const ValueType v3001, 
const ValueType v3010, const ValueType v0301, const ValueType v1030, const ValueType v0103,
const ValueType v2200, const ValueType v0220, const ValueType v0022, 
const ValueType v2002, const ValueType v2020, const ValueType v0202,
const ValueType v2110, const ValueType v0211, const ValueType v1021, const ValueType v1102, 
const ValueType v1210, const ValueType v0121, const ValueType v1012, const ValueType v2101, 
const ValueType v1120, const ValueType v0112, const ValueType v2011, const ValueType v1201,
const ValueType v1111)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);
  sA=(ValueType)4.0*sA;
  sB=(ValueType)4.0*sB;
  sC=(ValueType)4.0*sC;
  sD=(ValueType)4.0*sD;  

  if (sA>(ValueType)3.0) //1
  {
    result=(sA-(ValueType)3.0)*v4000+sB*v3100+sC*v3010+sD*v3001;
  }
  else if (sB>(ValueType)3.0) //2
  {
    result=sA*v1300+(sB-(ValueType)3.0)*v0400+sC*v0310+sD*v0301;
  }
  else if (sC>(ValueType)3.0) //3
  {
    result=sA*v1030+sB*v0130+(sC-(ValueType)3.0)*v0040+sD*v0031;
  }
  else if (sD>(ValueType)3.0) //4
  {
    result=sA*v1003+sB*v0103+sC*v0013+(sD-(ValueType)3.0)*v0004;
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)) //9
  {
    result=(sA-(ValueType)2.0)*v3100+(sB-(ValueType)1.0)*v2200+sC*v2110+sD*v2101;
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //10
  {
    result=sA*v1210+(sB-(ValueType)2.0)*v0310+(sC-(ValueType)1.0)*v0220+sD*v0211;
  }
  else if ((sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //11
  {
    result=sA*v1021+sB*v0121+(sC-(ValueType)2.0)*v0031+(sD-(ValueType)1.0)*v0022;
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //12
  {
    result=(sA-(ValueType)1.0)*v2002+sB*v1102+sC*v1012+(sD-(ValueType)2.0)*v1003;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)) //13
  {
    result=(sA-(ValueType)1.0)*v2200+(sB-(ValueType)2.0)*v1300+sC*v1210+sD*v1201;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //14
  {
    result=sA*v1120+(sB-(ValueType)1.0)*v0220+(sC-(ValueType)2.0)*v0130+sD*v0121;
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //15
  {
    result=sA*v1012+sB*v0112+(sC-(ValueType)1.0)*v0022+(sD-(ValueType)2.0)*v0013;
  }
  else if ((sA>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //16
  {
    result=(sA-(ValueType)2.0)*v3001+sB*v2101+sC*v2011+(sD-(ValueType)1.0)*v2002;
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //17
  {
    result=(sA-(ValueType)2.0)*v3010+sB*v2110+(sC-(ValueType)1.0)*v2020+sD*v2011;
  }
  else if ((sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //18
  {
    result=sA*v1201+(sB-(ValueType)2.0)*v0301+sC*v0211+(sD-(ValueType)1.0)*v0202;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //19
  {
    result=(sA-(ValueType)1.0)*v2020+sB*v1120+(sC-(ValueType)2.0)*v1030+sD*v1021;
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //20
  {
    result=sA*v1102+(sB-(ValueType)1.0)*v0202+sC*v0112+(sD-(ValueType)2.0)*v0103;
  }

  else if (sA>(ValueType)2.0) //25
  {
    if(axis=='0') //25.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sA)*v2110+sD*v3001+s1*v3100+s2*v3010; //25.1.1
        }
        else
        {
          result=sC*v2110+((ValueType)1.0-sB)*v3001+s1*v3100-s2*v2101; //25.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2110+((ValueType)1.0-sC)*v3001-s1*v2011+s2*v3010; //25.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2110+(sA-(ValueType)2.0)*v3001-s1*v2011-s2*v2101; //25.1.2
        };
      };
    }
    else if(axis=='1') //25.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v3010+sD*v2101+s1*v3100+s2*v2110; //25.2.1
        }
        else
        {
          result=sC*v3010+((ValueType)3.0-sA)*v2101+s1*v3100-s2*v3001; //25.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)2.0)*v3010+((ValueType)1.0-sC)*v2101-s1*v2011+s2*v2110; //25.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v3010+sB*v2101-s1*v2011-s2*v3001; //25.2.4
        };
      };
    }
    else //25.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v3100+sD*v2011+s1*v2110+s2*v3010; //25.3.1
        }
        else
        {
          result=(sA-(ValueType)2.0)*v3100+((ValueType)1.0-sB)*v2011+s1*v2110-s2*v2101; //25.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v3100+((ValueType)3.0-sA)*v2011-s1*v3001+s2*v3010; //25.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v3100+sC*v2011-s1*v3001-s2*v2101; //25.3.3
        };
      };
    };        
  }
  else if (sB>(ValueType)2.0) //26
  {
    if(axis=='0') //26.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0310+sD*v1201+s1*v1300+s2*v1210; //26.1.1
        }
        else
        {
          result=sC*v0310+((ValueType)3.0-sB)*v1201+s1*v1300-s2*v0301; //26.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v0310+((ValueType)1.0-sC)*v1201-s1*v0211+s2*v1210; //26.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0310+sA*v1201-s1*v0211-s2*v0301; //26.1.2
        };
      };
    }
    else if(axis=='1') //26.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sB)*v1210+sD*v0301+s1*v1300+s2*v0310; //26.2.1
        }
        else
        {
          result=sC*v1210+((ValueType)1.0-sA)*v0301+s1*v1300-s2*v1201; //26.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1210+((ValueType)1.0-sC)*v0301-s1*v0211+s2*v0310; //26.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1210+(sB-(ValueType)2.0)*v0301-s1*v0211-s2*v1201; //26.2.4
        };
      };
    }
    else //26.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1300+sD*v0211+s1*v0310+s2*v1210; //26.3.1
        }
        else
        {
          result=sA*v1300+((ValueType)3.0-sB)*v0211+s1*v0310-s2*v0301; //26.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v1300+((ValueType)1.0-sA)*v0211-s1*v1201+s2*v1210; //26.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1300+sC*v0211-s1*v1201-s2*v0301; //26.3.3
        };
      };
    };        
  }
  else if (sC>(ValueType)2.0) //27
  {
    if(axis=='0') //27.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0130+sD*v1021+s1*v1120+s2*v1030; //27.1.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v0130+((ValueType)1.0-sB)*v1021+s1*v1120-s2*v0121; //27.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0130+((ValueType)3.0-sC)*v1021-s1*v0031+s2*v1030; //27.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0130+sA*v1021-s1*v0031-s2*v0121; //27.1.2
        };
      };
    }
    else if(axis=='1') //27.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1030+sD*v0121+s1*v1120+s2*v0130; //27.2.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v1030+((ValueType)1.0-sA)*v0121+s1*v1120-s2*v1021; //27.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1030+((ValueType)3.0-sC)*v0121-s1*v0031+s2*v0130; //27.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1030+sB*v0121-s1*v0031-s2*v1021; //27.2.4
        };
      };
    }
    else //27.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sC)*v1120+sD*v0031+s1*v0130+s2*v1030; //27.3.1
        }
        else
        {
          result=sA*v1120+((ValueType)1.0-sB)*v0031+s1*v0130-s2*v0121; //27.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1120+((ValueType)1.0-sA)*v0031-s1*v1021+s2*v1030; //27.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1120+(sC-(ValueType)2.0)*v0031-s1*v1021-s2*v0121; //27.3.3
        };
      };
    };        
  }
  else if (sD>(ValueType)2.0) //28
  {
    if(axis=='0') //28.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0112+(sD-(ValueType)2.0)*v1003+s1*v1102+s2*v1012; //28.1.1
        }
        else
        {
          result=sC*v0112+((ValueType)1.0-sB)*v1003+s1*v1102-s2*v0103; //28.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0112+((ValueType)1.0-sC)*v1003-s1*v0013+s2*v1012; //28.1.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v0112+sA*v1003-s1*v0013-s2*v0103; //28.1.2
        };
      };
    }
    else if(axis=='1') //28.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1012+(sD-(ValueType)2.0)*v0103+s1*v1102+s2*v0112; //28.2.1
        }
        else
        {
          result=sC*v1012+((ValueType)1.0-sA)*v0103+s1*v1102-s2*v1003; //28.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1012+((ValueType)1.0-sC)*v0103-s1*v0013+s2*v0112; //28.2.2
        }
        else
        {
          result=((ValueType)3.0-sD)*v1012+sB*v0103-s1*v0013-s2*v1003; //28.2.4
        };
      };
    }
    else //28.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1102+(sD-(ValueType)2.0)*v0013+s1*v0112+s2*v1012; //28.3.1
        }
        else
        {
          result=sA*v1102+((ValueType)1.0-sB)*v0013+s1*v0112-s2*v0103; //28.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1102+((ValueType)1.0-sA)*v0013-s1*v1003+s2*v1012; //28.3.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v1102+sC*v0013-s1*v1003-s2*v0103; //28.3.3
        };
      };
    };        
  }

  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //21
  {
    result=((ValueType)2.0-sA)*v1111+((ValueType)1.0-sB)*v2011+((ValueType)1.0-sC)*v2101+((ValueType)1.0-sD)*v2110;
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //22
  {
    result=((ValueType)1.0-sA)*v0211+((ValueType)2.0-sB)*v1111+((ValueType)1.0-sC)*v1201+((ValueType)1.0-sD)*v1210;
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //23
  {
    result=((ValueType)1.0-sA)*v0121+((ValueType)1.0-sB)*v1021+((ValueType)2.0-sC)*v1111+((ValueType)1.0-sD)*v1120;
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //24
  {
    result=((ValueType)1.0-sA)*v0112+((ValueType)1.0-sB)*v1012+((ValueType)1.0-sC)*v1102+((ValueType)2.0-sD)*v1111;
  }

  else if ((sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //29
  {
    if(axis=='0') //29.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1210+sD*v2101+s1*v2200+s2*v2110; //29.1.1
        }
        else
        {
          result=sC*v1210+((ValueType)2.0-sB)*v2101+s1*v2200-s2*v1201; //29.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1210+((ValueType)1.0-sC)*v2101-s1*v1111+s2*v2110; //29.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1210+(sA-(ValueType)1.0)*v2101-s1*v1111-s2*v1201; //29.1.2
        };
      };
    }
    else if(axis=='1') //29.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v2110+sD*v1201+s1*v2200+s2*v1210; //29.2.1
        }
        else
        {
          result=sC*v2110+((ValueType)2.0-sA)*v1201+s1*v2200-s2*v2101; //29.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2110+((ValueType)1.0-sC)*v1201-s1*v1111+s2*v1210; //29.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v2110+(sB-(ValueType)1.0)*v1201-s1*v1111-s2*v2101; //29.2.4
        };
      };
    }
    else //29.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v2200+sD*v1111+s1*v1210+s2*v2110; //29.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2200+((ValueType)2.0-sB)*v1111+s1*v1210-s2*v1201; //29.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v2200+((ValueType)2.0-sA)*v1111-s1*v2101+s2*v2110; //29.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2200+sC*v1111-s1*v2101-s2*v1201; //29.3.3
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sD<(ValueType)1.0)) //30
  {
    if(axis=='0') //30.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0220+sD*v1111+s1*v1210+s2*v1120; //30.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v0220+((ValueType)2.0-sB)*v1111+s1*v1210-s2*v0211; //30.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v0220+((ValueType)2.0-sC)*v1111-s1*v0121+s2*v1120; //30.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0220+sA*v1111-s1*v0121-s2*v0211; //30.1.2
        };
      };
    }
    else if(axis=='1') //30.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v1120+sD*v0211+s1*v1210+s2*v0220; //30.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1120+((ValueType)1.0-sA)*v0211+s1*v1210-s2*v1111; //30.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1120+((ValueType)2.0-sC)*v0211-s1*v0121+s2*v0220; //30.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1120+(sB-(ValueType)1.0)*v0211-s1*v0121-s2*v1111; //30.2.4
        };
      };
    }
    else //30.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v1210+sD*v0121+s1*v0220+s2*v1120; //30.3.1
        }
        else
        {
          result=sA*v1210+((ValueType)2.0-sB)*v0121+s1*v0220-s2*v0211; //30.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1210+((ValueType)1.0-sA)*v0121-s1*v1111+s2*v1120; //30.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1210+(sC-(ValueType)1.0)*v0121-s1*v1111-s2*v0211; //30.3.3
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //31
  {
    if(axis=='0') //31.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0121+(sD-(ValueType)1.0)*v1012+s1*v1111+s2*v1021; //31.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v0121+((ValueType)1.0-sB)*v1012+s1*v1111-s2*v0112; //31.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0121+((ValueType)2.0-sC)*v1012-s1*v0022+s2*v1021; //31.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v0121+sA*v1012-s1*v0022-s2*v0112; //31.1.2
        };
      };
    }
    else if(axis=='1') //31.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1021+(sD-(ValueType)1.0)*v0112+s1*v1111+s2*v0121; //31.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1021+((ValueType)1.0-sA)*v0112+s1*v1111-s2*v1012; //31.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1021+((ValueType)2.0-sC)*v0112-s1*v0022+s2*v0121; //31.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v1021+sB*v0112-s1*v0022-s2*v1012; //31.2.4
        };
      };
    }
    else //31.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v1111+(sD-(ValueType)1.0)*v0022+s1*v0121+s2*v1021; //31.3.1
        }
        else
        {
          result=sA*v1111+((ValueType)1.0-sB)*v0022+s1*v0121-s2*v0112; //31.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1111+((ValueType)1.0-sA)*v0022-s1*v1012+s2*v1021; //31.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1111+(sC-(ValueType)1.0)*v0022-s1*v1012-s2*v0112; //31.3.3
        };
      };
    };        
  }
  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //32
  {
    if(axis=='0') //32.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1111+(sD-(ValueType)1.0)*v2002+s1*v2101+s2*v2011; //32.1.1
        }
        else
        {
          result=sC*v1111+((ValueType)1.0-sB)*v2002+s1*v2101-s2*v1102; //32.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1111+((ValueType)1.0-sC)*v2002-s1*v1012+s2*v2011; //32.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1111+(sA-(ValueType)1.0)*v2002-s1*v1012-s2*v1102; //32.1.2
        };
      };
    }
    else if(axis=='1') //32.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v2011+(sD-(ValueType)1.0)*v1102+s1*v2101+s2*v1111; //32.2.1
        }
        else
        {
          result=sC*v2011+((ValueType)2.0-sA)*v1102+s1*v2101-s2*v2002; //32.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2011+((ValueType)1.0-sC)*v1102-s1*v1012+s2*v1111; //32.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v2011+sB*v1102-s1*v1012-s2*v2002; //32.2.4
        };
      };
    }
    else //32.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v2101+(sD-(ValueType)1.0)*v1012+s1*v1111+s2*v2011; //32.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2101+((ValueType)1.0-sB)*v1012+s1*v1111-s2*v1102; //32.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2101+((ValueType)2.0-sA)*v1012-s1*v2002+s2*v2011; //32.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v2101+sC*v1012-s1*v2002-s2*v1102; //32.3.3
        };
      };
    };        
  }
  else if ((sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //33
  {
    if(axis=='0') //33.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1120+sD*v2011+s1*v2200+s2*v2110; //33.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1120+((ValueType)1.0-sB)*v2011+s1*v2200-s2*v1201; //33.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1120+((ValueType)2.0-sC)*v2011-s1*v1111+s2*v2110; //33.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1120+(sA-(ValueType)1.0)*v2011-s1*v1111-s2*v1201; //33.1.2
        };
      };
    }
    else if(axis=='1') //33.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v2020+sD*v1111+s1*v2110+s2*v1120; //33.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v2020+((ValueType)2.0-sA)*v1111+s1*v2110-s2*v2011; //33.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2020+((ValueType)2.0-sC)*v1111-s1*v1021+s2*v1120; //33.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v2020+sB*v1111-s1*v1021-s2*v2011; //33.2.4
        };
      };
    }
    else //33.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v2110+sD*v1021+s1*v1120+s2*v2020; //33.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2110+((ValueType)1.0-sB)*v1021+s1*v1120-s2*v1111; //33.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2110+((ValueType)2.0-sA)*v1021-s1*v2011+s2*v2020; //33.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2110+(sC-(ValueType)1.0)*v1021-s1*v2011-s2*v1111; //33.3.3
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //34
  {
    if(axis=='0') //34.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0211+(sD-(ValueType)1.0)*v1102+s1*v1201+s2*v1111; //34.1.1
        }
        else
        {
          result=sC*v0211+((ValueType)2.0-sB)*v1102+s1*v1201-s2*v0202; //34.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v0211+((ValueType)1.0-sC)*v1102-s1*v0112+s2*v1111; //34.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v0211+sA*v1102-s1*v0112-s2*v0202; //34.1.2
        };
      };
    }
    else if(axis=='1') //34.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v1111+(sD-(ValueType)1.0)*v0202+s1*v1201+s2*v0211; //34.2.1
        }
        else
        {
          result=sC*v1111+((ValueType)1.0-sA)*v0202+s1*v1201-s2*v1102; //34.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1111+((ValueType)1.0-sC)*v0202-s1*v0112+s2*v0211; //34.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v1111+(sB-(ValueType)1.0)*v0202-s1*v0112-s2*v1102; //34.2.4
        };
      };
    }
    else //34.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1201+(sD-(ValueType)1.0)*v0112+s1*v0211+s2*v1111; //34.3.1
        }
        else
        {
          result=sA*v1201+((ValueType)2.0-sB)*v0112+s1*v0211-s2*v0202; //34.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1201+((ValueType)1.0-sA)*v0112-s1*v1102+s2*v1111; //34.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1201+sC*v0112-s1*v1102-s2*v0202; //34.3.3
        };
      };
    };        
  }

  else if (sA<(ValueType)1.0) //5
  {
    result=sA*v1111+(sB-(ValueType)1.0)*v0211+(sC-(ValueType)1.0)*v0121+(sD-(ValueType)1.0)*v0112;
  }
  else if (sB<(ValueType)1.0) //6
  {
    result=(sA-(ValueType)1.0)*v2011+sB*v1111+(sC-(ValueType)1.0)*v1021+(sD-(ValueType)1.0)*v1012;
  }
  else if (sC<(ValueType)1.0) //7
  {
    result=(sA-(ValueType)1.0)*v2101+(sB-(ValueType)1.0)*v1201+sC*v1111+(sD-(ValueType)1.0)*v1102;
  }
  else //8
  {
    result=(sA-(ValueType)1.0)*v2110+(sB-(ValueType)1.0)*v1210+(sC-(ValueType)1.0)*v1120+sD*v1111;
  };

  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Mono_4(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v4000, const ValueType v0400, const ValueType v0040, const ValueType v0004,
const ValueType v3100, const ValueType v0310, const ValueType v0031, const ValueType v1003, 
const ValueType v1300, const ValueType v0130, const ValueType v0013, const ValueType v3001, 
const ValueType v3010, const ValueType v0301, const ValueType v1030, const ValueType v0103,
const ValueType v2200, const ValueType v0220, const ValueType v0022, 
const ValueType v2002, const ValueType v2020, const ValueType v0202,
const ValueType v2110, const ValueType v0211, const ValueType v1021, const ValueType v1102, 
const ValueType v1210, const ValueType v0121, const ValueType v1012, const ValueType v2101, 
const ValueType v1120, const ValueType v0112, const ValueType v2011, const ValueType v1201,
const ValueType v1111)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;
  ValueType Vmin, Vmax;
  ValueType w4000, w0400, w0040, w0004;
  ValueType w3100, w0310, w0031, w1003; 
  ValueType w1300, w0130, w0013, w3001; 
  ValueType w3010, w0301, w1030, w0103;
  ValueType w2200, w0220, w0022, w2002, w2020, w0202;
  ValueType w2110, w0211, w1021, w1102; 
  ValueType w1210, w0121, w1012, w2101; 
  ValueType w1120, w0112, w2011, w1201;
  ValueType w1111;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  w4000=sA*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)3.0)/(ValueType)3.0;
  w0400=sB*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)3.0)/(ValueType)3.0;
  w0040=sC*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)3.0)/(ValueType)3.0;
  w0004=sD*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*((ValueType)4.0*sD-(ValueType)3.0)/(ValueType)3.0;

  w3100=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0310=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0031=sC*sD*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w1003=sD*sA*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w1300=sB*sA*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0130=sC*sB*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0013=sD*sC*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w3001=sA*sD*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w3010=sA*sC*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)2.0*sA-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0301=sB*sD*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)2.0*sB-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w1030=sC*sA*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)2.0*sC-(ValueType)1.0)*(ValueType)(16.0/3.0);
  w0103=sD*sB*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)2.0*sD-(ValueType)1.0)*(ValueType)(16.0/3.0);

  w2200=sA*sB*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)4.0;
  w0220=sB*sC*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)4.0;
  w0022=sC*sD*((ValueType)4.0*sC-(ValueType)1.0)*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)4.0;
  w2002=sD*sA*((ValueType)4.0*sD-(ValueType)1.0)*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)4.0;
  w2020=sA*sC*((ValueType)4.0*sA-(ValueType)1.0)*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)4.0;
  w0202=sB*sD*((ValueType)4.0*sB-(ValueType)1.0)*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)4.0;

  w2110=sA*sB*sC*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w0211=sB*sC*sD*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;
  w1021=sC*sD*sA*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;
  w1102=sD*sA*sB*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)32.0;
  w1210=sB*sC*sA*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;
  w0121=sC*sD*sB*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;
  w1012=sD*sA*sC*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)32.0;
  w2101=sA*sB*sD*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w1120=sC*sA*sB*((ValueType)4.0*sC-(ValueType)1.0)*(ValueType)32.0;
  w0112=sD*sB*sC*((ValueType)4.0*sD-(ValueType)1.0)*(ValueType)32.0;
  w2011=sA*sC*sD*((ValueType)4.0*sA-(ValueType)1.0)*(ValueType)32.0;
  w1201=sB*sD*sA*((ValueType)4.0*sB-(ValueType)1.0)*(ValueType)32.0;

  w1111=sA*sB*sC*sD*(ValueType)256.0;

  result=v4000*w4000+v0400*w0400+v0040*w0040+v0004*w0004+
         v3100*w3100+v0310*w0310+v0031*w0031+v1003*w1003+
         v1300*w1300+v0130*w0130+v0013*w0013+v3001*w3001+
         v3010*w3010+v0301*w0301+v1030*w1030+v0103*w0103+
         v2200*w2200+v0220*w0220+v0022*w0022+v2002*w2002+v2020*w2020+v0202*w0202+
         v2110*w2110+v0211*w0211+v1021*w1021+v1102*w1102+
         v1210*w1210+v0121*w0121+v1012*w1012+v2101*w2101+
         v1120*w1120+v0112*w0112+v2011*w2011+v1201*w1201+
         v1111*w1111;

  sA=(ValueType)4.0*sA;
  sB=(ValueType)4.0*sB;
  sC=(ValueType)4.0*sC;
  sD=(ValueType)4.0*sD;  

  if (sA>(ValueType)3.0) //1
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4000, v3100, v3010, v3001);
  }
  else if (sB>(ValueType)3.0) //2
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1300, v0400, v0310, v0301);
  }
  else if (sC>(ValueType)3.0) //3
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1030, v0130, v0040, v0031);
  }
  else if (sD>(ValueType)3.0) //4
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1003, v0103, v0013, v0004);
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)) //9
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3100, v2200, v2110, v2101);
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //10
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0310, v0220, v0211);
  }
  else if ((sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //11
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1021, v0121, v0031, v0022);
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //12
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2002, v1102, v1012, v1003);
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)) //13
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2200, v1300, v1210, v1201);
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //14
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0220, v0130, v0121);
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //15
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1012, v0112, v0022, v0013);
  }
  else if ((sA>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //16
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3001, v2101, v2011, v2002);
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //17
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3010, v2110, v2020, v2011);
  }
  else if ((sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //18
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1201, v0301, v0211, v0202);
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //19
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2020, v1120, v1030, v1021);
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //20
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1102, v0202, v0112, v0103);
  }

  else if (sA>(ValueType)2.0) //25
  {
    if(axis=='0') //25.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v3001, v3100, v3010); //25.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v3001, v3100, v2101); //25.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v3001, v2011, v3010); //25.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v3001, v2011, v2101); //25.1.2
        };
      };
    }
    else if(axis=='1') //25.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3010, v2101, v3100, v2110); //25.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3010, v2101, v3100, v3001); //25.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3010, v2101, v2011, v2110); //25.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3010, v2101, v2011, v3001); //25.2.4
        };
      };
    }
    else //25.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3100, v2011, v2110, v3010); //25.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3100, v2011, v2110, v2101); //25.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3100, v2011, v3001, v3010); //25.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3100, v2011, v3001, v2101); //25.3.3
        };
      };
    };        
  }
  else if (sB>(ValueType)2.0) //26
  {
    if(axis=='0') //26.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0310, v1201, v1300, v1210); //26.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0310, v1201, v1300, v0301); //26.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0310, v1201, v0211, v1210); //26.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0310, v1201, v0211, v0301); //26.1.2
        };
      };
    }
    else if(axis=='1') //26.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0301, v1300, v0310); //26.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0301, v1300, v1201); //26.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0301, v0211, v0310); //26.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0301, v0211, v1201); //26.2.4
        };
      };
    }
    else //26.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1300, v0211, v0310, v1210); //26.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1300, v0211, v0310, v0301); //26.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1300, v0211, v1201, v1210); //26.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1300, v0211, v1201, v0301); //26.3.3
        };
      };
    };        
  }
  else if (sC>(ValueType)2.0) //27
  {
    if(axis=='0') //27.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0130, v1021, v1120, v1030); //27.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0130, v1021, v1120, v0121); //27.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0130, v1021, v0031, v1030); //27.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0130, v1021, v0031, v0121); //27.1.2
        };
      };
    }
    else if(axis=='1') //27.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1030, v0121, v1120, v0130); //27.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1030, v0121, v1120, v1021); //27.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1030, v0121, v0031, v0130); //27.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1030, v0121, v0031, v1021); //27.2.4
        };
      };
    }
    else //27.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0031, v0130, v1030); //27.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0031, v0130, v0121); //27.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0031, v1021, v1030); //27.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0031, v1021, v0121); //27.3.3
        };
      };
    };        
  }
  else if (sD>(ValueType)2.0) //28
  {
    if(axis=='0') //28.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0112, v1003, v1102, v1012); //28.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0112, v1003, v1102, v0103); //28.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0112, v1003, v0013, v1012); //28.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0112, v1003, v0013, v0103); //28.1.2
        };
      };
    }
    else if(axis=='1') //28.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1012, v0103, v1102, v0112); //28.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1012, v0103, v1102, v1003); //28.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1012, v0103, v0013, v0112); //28.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1012, v0103, v0013, v1003); //28.2.4
        };
      };
    }
    else //28.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1102, v0013, v0112, v1012); //28.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1102, v0013, v0112, v0103); //28.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1102, v0013, v1003, v1012); //28.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1102, v0013, v1003, v0103); //28.3.3
        };
      };
    };        
  }

  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //21
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v2011, v2101, v2110);
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //22
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0211, v1111, v1201, v1210);
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //23
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0121, v1021, v1111, v1120);
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //24
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0112, v1012, v1102, v1111);
  }

  else if ((sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //29
  {
    if(axis=='0') //29.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v2101, v2200, v2110); //29.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v2101, v2200, v1201); //29.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v2101, v1111, v2110); //29.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v2101, v1111, v1201); //29.1.2
        };
      };
    }
    else if(axis=='1') //29.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1201, v2200, v1210); //29.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1201, v2200, v2101); //29.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1201, v1111, v1210); //29.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1201, v1111, v2101); //29.2.4
        };
      };
    }
    else //29.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2200, v1111, v1210, v2110); //29.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2200, v1111, v1210, v1201); //29.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2200, v1111, v2101, v2110); //29.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2200, v1111, v2101, v1201); //29.3.3
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sD<(ValueType)1.0)) //30
  {
    if(axis=='0') //30.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0220, v1111, v1210, v1120); //30.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0220, v1111, v1210, v0211); //30.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0220, v1111, v0121, v1120); //30.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0220, v1111, v0121, v0211); //30.1.2
        };
      };
    }
    else if(axis=='1') //30.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0211, v1210, v0220); //30.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0211, v1210, v1111); //30.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0211, v0121, v0220); //30.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v0211, v0121, v1111); //30.2.4
        };
      };
    }
    else //30.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0121, v0220, v1120); //30.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0121, v0220, v0211); //30.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0121, v1111, v1120); //30.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1210, v0121, v1111, v0211); //30.3.3
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //31
  {
    if(axis=='0') //31.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0121, v1012, v1111, v1021); //31.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0121, v1012, v1111, v0112); //31.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0121, v1012, v0022, v1021); //31.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0121, v1012, v0022, v0112); //31.1.2
        };
      };
    }
    else if(axis=='1') //31.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1021, v0112, v1111, v0121); //31.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1021, v0112, v1111, v1012); //31.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1021, v0112, v0022, v0121); //31.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1021, v0112, v0022, v1012); //31.2.4
        };
      };
    }
    else //31.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0022, v0121, v1021); //31.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0022, v0121, v0112); //31.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0022, v1012, v1021); //31.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0022, v1012, v0112); //31.3.3
        };
      };
    };        
  }
  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //32
  {
    if(axis=='0') //32.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v2002, v2101, v2011); //32.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v2002, v2101, v1102); //32.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v2002, v1012, v2011); //32.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v2002, v1012, v1102); //32.1.2
        };
      };
    }
    else if(axis=='1') //32.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2011, v1102, v2101, v1111); //32.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2011, v1102, v2101, v2002); //32.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2011, v1102, v1012, v1111); //32.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2011, v1102, v1012, v2002); //32.2.4
        };
      };
    }
    else //32.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2101, v1012, v1111, v2011); //32.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2101, v1012, v1111, v1102); //32.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2101, v1012, v2002, v2011); //32.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2101, v1012, v2002, v1102); //32.3.3
        };
      };
    };        
  }
  else if ((sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //33
  {
    if(axis=='0') //33.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v2011, v2200, v2110); //33.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v2011, v2200, v1201); //33.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v2011, v1111, v2110); //33.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1120, v2011, v1111, v1201); //33.1.2
        };
      };
    }
    else if(axis=='1') //33.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2020, v1111, v2110, v1120); //33.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2020, v1111, v2110, v2011); //33.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2020, v1111, v1021, v1120); //33.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2020, v1111, v1021, v2011); //33.2.4
        };
      };
    }
    else //33.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1021, v1120, v2020); //33.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1021, v1120, v1111); //33.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1021, v2011, v2020); //33.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1021, v2011, v1111); //33.3.3
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //34
  {
    if(axis=='0') //34.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0211, v1102, v1201, v1111); //34.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0211, v1102, v1201, v0202); //34.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0211, v1102, v0112, v1111); //34.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0211, v1102, v0112, v0202); //34.1.2
        };
      };
    }
    else if(axis=='1') //34.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0202, v1201, v0211); //34.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0202, v1201, v1102); //34.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0202, v0112, v0211); //34.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0202, v0112, v1102); //34.2.4
        };
      };
    }
    else //34.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1201, v0112, v0211, v1111); //34.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1201, v0112, v0211, v0202); //34.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1201, v0112, v1102, v1111); //34.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1201, v0112, v1102, v0202); //34.3.3
        };
      };
    };        
  }

  else if (sA<(ValueType)1.0) //5
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1111, v0211, v0121, v0112);
  }
  else if (sB<(ValueType)1.0) //6
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2011, v1111, v1021, v1012);
  }
  else if (sC<(ValueType)1.0) //7
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2101, v1201, v1111, v1102);
  }
  else //8
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2110, v1210, v1120, v1111);
  };

  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };

  return(result);
}

////////////////////////////////////////////////N=5//////////////////////////////////////////

//2D

template <typename ValueType>
ValueType Interpolate_2D_Power_5(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v500, const ValueType v050, const ValueType v005,
const ValueType v410, const ValueType v041, const ValueType v104, 
const ValueType v140, const ValueType v014, const ValueType v401,
const ValueType v320, const ValueType v032, const ValueType v203, 
const ValueType v230, const ValueType v023, const ValueType v302,
const ValueType v311, const ValueType v131, const ValueType v113,
const ValueType v221, const ValueType v122, const ValueType v212)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType w500, w050, w005;
  ValueType w410, w041, w104; 
  ValueType w140, w014, w401;
  ValueType w320, w032, w203; 
  ValueType w230, w023, w302;
  ValueType w311, w131, w113;
  ValueType w221, w122, w212;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);

  w500=sA*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*((ValueType)5.0*sA-(ValueType)4.0)/(ValueType)24.0;
  w050=sB*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*((ValueType)5.0*sB-(ValueType)4.0)/(ValueType)24.0;
  w005=sC*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*((ValueType)5.0*sC-(ValueType)4.0)/(ValueType)24.0;

  w410=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w041=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w104=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w140=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w014=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w401=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);

  w320=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w032=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w203=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w230=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w023=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w302=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w311=sA*sB*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w131=sB*sC*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w113=sC*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);

  w221=sA*sB*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w122=sB*sC*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w212=sC*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);

  result=v500*w500+v050*w050+v005*w005+
         v410*w410+v041*w041+v104*w104+
         v140*w140+v014*w014+v401*w401+
         v320*w320+v032*w032+v203*w203+
         v230*w230+v023*w023+v302*w302+
         v311*w311+v131*w131+v113*w113+
         v221*w221+v122*w122+v212*w212;
  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Linear_5(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v500, const ValueType v050, const ValueType v005,
const ValueType v410, const ValueType v041, const ValueType v104, 
const ValueType v140, const ValueType v014, const ValueType v401,
const ValueType v320, const ValueType v032, const ValueType v203, 
const ValueType v230, const ValueType v023, const ValueType v302,
const ValueType v311, const ValueType v131, const ValueType v113,
const ValueType v221, const ValueType v122, const ValueType v212)
{
  ValueType result;
  ValueType sA, sB, sC;
  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);
  sA=(ValueType)5.0*sA;
  sB=(ValueType)5.0*sB;
  sC=(ValueType)5.0*sC;
  if (sA>(ValueType)4.0) //1
  {
    result=(sA-(ValueType)4.0)*v500+sB*v410+sC*v401;
  }
  else if (sB>(ValueType)4.0) //2
  {
    result=sA*v140+(sB-(ValueType)4.0)*v050+sC*v041;
  }
  else if (sC>(ValueType)4.0) //3
  {
    result=sA*v104+sB*v014+(sC-(ValueType)4.0)*v005;
  }

  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //16
  {
    result=((ValueType)4.0-sA)*v311+((ValueType)1.0-sB)*v401+((ValueType)1.0-sC)*v410;
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //17
  {
    result=((ValueType)1.0-sA)*v041+((ValueType)4.0-sB)*v131+((ValueType)1.0-sC)*v140;
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //18
  {
    result=((ValueType)1.0-sA)*v014+((ValueType)1.0-sB)*v104+((ValueType)4.0-sC)*v113;
  }

  else if ((sA>(ValueType)3.0)&&(sC<(ValueType)1.0)) //4
  {
    result=(sA-(ValueType)3.0)*v410+(sB-(ValueType)1.0)*v320+sC*v311;
  }
  else if ((sA<(ValueType)1.0)&&(sB>(ValueType)3.0)) //5
  {
    result=sA*v131+(sB-(ValueType)3.0)*v041+(sC-(ValueType)1.0)*v032;
  }
  else if ((sB<(ValueType)1.0)&&(sC>(ValueType)3.0)) //6
  {
    result=(sA-(ValueType)1.0)*v203+sB*v113+(sC-(ValueType)3.0)*v104;
  }
  else if ((sB>(ValueType)3.0)&&(sC<(ValueType)1.0)) //7
  {
    result=(sA-(ValueType)1.0)*v230+(sB-(ValueType)3.0)*v140+sC*v131;
  }
  else if ((sA<(ValueType)1.0)&&(sC>(ValueType)3.0)) //8
  {
    result=sA*v113+(sB-(ValueType)1.0)*v023+(sC-(ValueType)3.0)*v014;
  }
  else if ((sA>(ValueType)3.0)&&(sB<(ValueType)1.0)) //9
  {
    result=(sA-(ValueType)3.0)*v401+sB*v311+(sC-(ValueType)1.0)*v302;
  }
  else if ((sA>(ValueType)2.0)&&(sB>(ValueType)2.0)) //10
  {
    result=(sA-(ValueType)2.0)*v320+(sB-(ValueType)2.0)*v230+sC*v221;
  }
  else if ((sB>(ValueType)2.0)&&(sC>(ValueType)2.0)) //11
  {
    result=sA*v122+(sB-(ValueType)2.0)*v032+(sC-(ValueType)2.0)*v023;
  }
  else if ((sA>(ValueType)2.0)&&(sC>(ValueType)2.0)) //12
  {
    result=(sA-(ValueType)2.0)*v302+sB*v212+(sC-(ValueType)2.0)*v203;
  }

  else if ((sA>(ValueType)2.0)&&(sC<(ValueType)1.0)) //19
  {
    result=((ValueType)3.0-sA)*v221+((ValueType)2.0-sB)*v311+((ValueType)1.0-sC)*v320;
  }
  else if ((sA<(ValueType)1.0)&&(sB>(ValueType)2.0)) //20
  {
    result=((ValueType)1.0-sA)*v032+((ValueType)3.0-sB)*v122+((ValueType)2.0-sC)*v131;
  }
  else if ((sB<(ValueType)1.0)&&(sC>(ValueType)2.0)) //21
  {
    result=((ValueType)2.0-sA)*v113+((ValueType)1.0-sB)*v203+((ValueType)3.0-sC)*v212;
  }
  else if ((sA<(ValueType)1.0)&&(sC>(ValueType)2.0)) //22
  {
    result=((ValueType)1.0-sA)*v023+((ValueType)2.0-sB)*v113+((ValueType)3.0-sC)*v122;
  }
  else if ((sA>(ValueType)2.0)&&(sB<(ValueType)1.0)) //23
  {
    result=((ValueType)3.0-sA)*v212+((ValueType)1.0-sB)*v302+((ValueType)2.0-sC)*v311;
  }
  else if ((sB>(ValueType)2.0)&&(sC<(ValueType)1.0)) //24
  {
    result=((ValueType)2.0-sA)*v131+((ValueType)3.0-sB)*v221+((ValueType)1.0-sC)*v230;
  }

  else if (sA>(ValueType)2.0) //13
  {
    result=(sA-(ValueType)2.0)*v311+(sB-(ValueType)1.0)*v221+(sC-(ValueType)1.0)*v212;
  }
  else if (sB>(ValueType)2.0) //14
  {
    result=(sA-(ValueType)1.0)*v221+(sB-(ValueType)2.0)*v131+(sC-(ValueType)1.0)*v122;
  }
  else if (sC>(ValueType)2.0) //15
  {
    result=(sA-(ValueType)1.0)*v212+(sB-(ValueType)1.0)*v122+(sC-(ValueType)2.0)*v113;
  }

  else //25
  {
    result=((ValueType)2.0-sA)*v122+((ValueType)2.0-sB)*v212+((ValueType)2.0-sC)*v221;
  };

  return(result);
}

template <typename ValueType>
ValueType Interpolate_2D_Mono_5(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v500, const ValueType v050, const ValueType v005,
const ValueType v410, const ValueType v041, const ValueType v104, 
const ValueType v140, const ValueType v014, const ValueType v401,
const ValueType v320, const ValueType v032, const ValueType v203, 
const ValueType v230, const ValueType v023, const ValueType v302,
const ValueType v311, const ValueType v131, const ValueType v113,
const ValueType v221, const ValueType v122, const ValueType v212)
{
  ValueType result;
  ValueType sA, sB, sC;
  ValueType Vmin, Vmax;
  ValueType w500, w050, w005;
  ValueType w410, w041, w104; 
  ValueType w140, w014, w401;
  ValueType w320, w032, w203; 
  ValueType w230, w023, w302;
  ValueType w311, w131, w113;
  ValueType w221, w122, w212;

  interpolationInitial_2D<ValueType>(x, y, 
                                     xA, yA, 
                                     xB, yB, 
                                     xC, yC, 
                                     &sA, &sB, &sC);

  w500=sA*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*((ValueType)5.0*sA-(ValueType)4.0)/(ValueType)24.0;
  w050=sB*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*((ValueType)5.0*sB-(ValueType)4.0)/(ValueType)24.0;
  w005=sC*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*((ValueType)5.0*sC-(ValueType)4.0)/(ValueType)24.0;

  w410=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w041=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w104=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w140=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w014=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w401=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);

  w320=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w032=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w203=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w230=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w023=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w302=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w311=sA*sB*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w131=sB*sC*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w113=sC*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);

  w221=sA*sB*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w122=sB*sC*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w212=sC*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);

  result=v500*w500+v050*w050+v005*w005+
         v410*w410+v041*w041+v104*w104+
         v140*w140+v014*w014+v401*w401+
         v320*w320+v032*w032+v203*w203+
         v230*w230+v023*w023+v302*w302+
         v311*w311+v131*w131+v113*w113+
         v221*w221+v122*w122+v212*w212;

  sA=(ValueType)5.0*sA;
  sB=(ValueType)5.0*sB;
  sC=(ValueType)5.0*sC;

  if (sA>(ValueType)4.0) //1
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v500, v410, v401);
  }
  else if (sB>(ValueType)4.0) //2
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v140, v050, v041);
  }
  else if (sC>(ValueType)4.0) //3
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v104, v014, v005);
  }

  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //16
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v311, v401, v410);
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //17
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v041, v131, v140);
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //18
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v014, v104, v113);
  }

  else if ((sA>(ValueType)3.0)&&(sC<(ValueType)1.0)) //4
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v410, v320, v311);
  }
  else if ((sA<(ValueType)1.0)&&(sB>(ValueType)3.0)) //5
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v131, v041, v032);
  }
  else if ((sB<(ValueType)1.0)&&(sC>(ValueType)3.0)) //6
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v203, v113, v104);
  }
  else if ((sB>(ValueType)3.0)&&(sC<(ValueType)1.0)) //7
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v230, v140, v131);
  }
  else if ((sA<(ValueType)1.0)&&(sC>(ValueType)3.0)) //8
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v113, v023, v014);
  }
  else if ((sA>(ValueType)3.0)&&(sB<(ValueType)1.0)) //9
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v401, v311, v302);
  }
  else if ((sA>(ValueType)2.0)&&(sB>(ValueType)2.0)) //10
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v320, v230, v221);
  }
  else if ((sB>(ValueType)2.0)&&(sC>(ValueType)2.0)) //11
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v122, v032, v023);
  }
  else if ((sA>(ValueType)2.0)&&(sC>(ValueType)2.0)) //12
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v302, v212, v203);
  }

  else if ((sA>(ValueType)2.0)&&(sC<(ValueType)1.0)) //19
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v221, v311, v320);
  }
  else if ((sA<(ValueType)1.0)&&(sB>(ValueType)2.0)) //20
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v032, v122, v131);
  }
  else if ((sB<(ValueType)1.0)&&(sC>(ValueType)2.0)) //21
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v113, v203, v212);
  }
  else if ((sA<(ValueType)1.0)&&(sC>(ValueType)2.0)) //22
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v023, v113, v122);
  }
  else if ((sA>(ValueType)2.0)&&(sB<(ValueType)1.0)) //23
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v212, v302, v311);
  }
  else if ((sB>(ValueType)2.0)&&(sC<(ValueType)1.0)) //24
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v131, v221, v230);
  }

  else if (sA>(ValueType)2.0) //13
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v311, v221, v212);
  }
  else if (sB>(ValueType)2.0) //14
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v221, v131, v122);
  }
  else if (sC>(ValueType)2.0) //15
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v212, v122, v113);
  }

  else //25
  {
    interpolationMinMax_2D<ValueType>(&Vmin, &Vmax, v122, v212, v221);
  };

  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };
  return(result);
}

//3D

template <typename ValueType>
ValueType Interpolate_3D_Power_5(
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v5000, const ValueType v0500, const ValueType v0050, const ValueType v0005,
const ValueType v4100, const ValueType v0410, const ValueType v0041, const ValueType v1004, 
const ValueType v1400, const ValueType v0140, const ValueType v0014, const ValueType v4001, 
const ValueType v4010, const ValueType v0401, const ValueType v1040, const ValueType v0104,
const ValueType v3200, const ValueType v0320, const ValueType v0032, const ValueType v2003, 
const ValueType v2300, const ValueType v0230, const ValueType v0023, const ValueType v3002, 
const ValueType v3020, const ValueType v0302, const ValueType v2030, const ValueType v0203,
const ValueType v3110, const ValueType v0311, const ValueType v1031, const ValueType v1103, 
const ValueType v1310, const ValueType v0131, const ValueType v1013, const ValueType v3101, 
const ValueType v1130, const ValueType v0113, const ValueType v3011, const ValueType v1301,
const ValueType v2210, const ValueType v0221, const ValueType v1022, const ValueType v2102,
const ValueType v1220, const ValueType v0122, const ValueType v2012, const ValueType v2201, 
const ValueType v2120, const ValueType v0212, const ValueType v2021, const ValueType v1202, 
const ValueType v2111, const ValueType v1211, const ValueType v1121, const ValueType v1112)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType w5000, w0500, w0050, w0005;
  ValueType w4100, w0410, w0041, w1004; 
  ValueType w1400, w0140, w0014, w4001; 
  ValueType w4010, w0401, w1040, w0104;
  ValueType w3200, w0320, w0032, w2003; 
  ValueType w2300, w0230, w0023, w3002; 
  ValueType w3020, w0302, w2030, w0203;
  ValueType w3110, w0311, w1031, w1103; 
  ValueType w1310, w0131, w1013, w3101; 
  ValueType w1130, w0113, w3011, w1301;
  ValueType w2210, w0221, w1022, w2102;
  ValueType w1220, w0122, w2012, w2201; 
  ValueType w2120, w0212, w2021, w1202; 
  ValueType w2111, w1211, w1121, w1112;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  w5000=sA*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*((ValueType)5.0*sA-(ValueType)4.0)/(ValueType)24.0;
  w0500=sB*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*((ValueType)5.0*sB-(ValueType)4.0)/(ValueType)24.0;
  w0050=sC*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*((ValueType)5.0*sC-(ValueType)4.0)/(ValueType)24.0;
  w0005=sD*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*((ValueType)5.0*sD-(ValueType)4.0)/(ValueType)24.0;

  w4100=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0410=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0041=sC*sD*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w1004=sD*sA*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w1400=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0140=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0014=sD*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w4001=sA*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w4010=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0401=sB*sD*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w1040=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0104=sD*sB*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*(ValueType)(25.0/24.0);

  w3200=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0320=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0032=sC*sD*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w2003=sD*sA*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w2300=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0230=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0023=sD*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w3002=sA*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w3020=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0302=sB*sD*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w2030=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0203=sD*sB*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w3110=sA*sB*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w0311=sB*sC*sD*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1031=sC*sD*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1103=sD*sA*sB*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1310=sB*sC*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w0131=sC*sD*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1013=sD*sA*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w3101=sA*sB*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1130=sC*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w0113=sD*sB*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w3011=sA*sC*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1301=sB*sD*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);

  w2210=sA*sB*sC*(5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w0221=sB*sC*sD*(5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w1022=sC*sD*sA*(5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2102=sD*sA*sB*(5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w1220=sB*sC*sA*(5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w0122=sC*sD*sB*(5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2012=sD*sA*sC*(5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2201=sA*sB*sD*(5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2120=sC*sA*sB*(5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w0212=sD*sB*sC*(5.0*sD-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2021=sA*sC*sD*(5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w1202=sB*sD*sA*(5.0*sB-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(125.0/4.0);

  w2111=sA*sB*sC*sD*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(625.0/2.0);
  w1211=sB*sC*sD*sA*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(625.0/2.0);
  w1121=sC*sD*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(625.0/2.0);
  w1112=sD*sA*sB*sC*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(625.0/2.0);

  result=v5000*w5000+v0500*w0500+v0050*w0050+v0005*w0005+
         v4100*w4100+v0410*w0410+v0041*w0041+v1004*w1004+
         v1400*w1400+v0140*w0140+v0014*w0014+v4001*w4001+
         v4010*w4010+v0401*w0401+v1040*w1040+v0104*w0104+
         v3200*w3200+v0320*w0320+v0032*w0032+v2003*w2003+
         v2300*w2300+v0230*w0230+v0023*w0023+v3002*w3002+
         v3020*w3020+v0302*w0302+v2030*w2030+v0203*w0203+
         v3110*w3110+v0311*w0311+v1031*w1031+v1103*w1103+
         v1310*w1310+v0131*w0131+v1013*w1013+v3101*w3101+
         v1130*w1130+v0113*w0113+v3011*w3011+v1301*w1301+
         v2210*w2210+v0221*w0221+v1022*w1022+v2102*w2102+
         v1220*w1220+v0122*w0122+v2012*w2012+v2201*w2201+
         v2120*w2120+v0212*w0212+v2021*w2021+v1202*w1202+
         v2111*w2111+v1211*w1211+v1121*w1121+v1112*w1112;
  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Linear_5(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v5000, const ValueType v0500, const ValueType v0050, const ValueType v0005,
const ValueType v4100, const ValueType v0410, const ValueType v0041, const ValueType v1004, 
const ValueType v1400, const ValueType v0140, const ValueType v0014, const ValueType v4001, 
const ValueType v4010, const ValueType v0401, const ValueType v1040, const ValueType v0104,
const ValueType v3200, const ValueType v0320, const ValueType v0032, const ValueType v2003, 
const ValueType v2300, const ValueType v0230, const ValueType v0023, const ValueType v3002, 
const ValueType v3020, const ValueType v0302, const ValueType v2030, const ValueType v0203,
const ValueType v3110, const ValueType v0311, const ValueType v1031, const ValueType v1103, 
const ValueType v1310, const ValueType v0131, const ValueType v1013, const ValueType v3101, 
const ValueType v1130, const ValueType v0113, const ValueType v3011, const ValueType v1301,
const ValueType v2210, const ValueType v0221, const ValueType v1022, const ValueType v2102,
const ValueType v1220, const ValueType v0122, const ValueType v2012, const ValueType v2201, 
const ValueType v2120, const ValueType v0212, const ValueType v2021, const ValueType v1202, 
const ValueType v2111, const ValueType v1211, const ValueType v1121, const ValueType v1112)
{



  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);
  sA=(ValueType)5.0*sA;
  sB=(ValueType)5.0*sB;
  sC=(ValueType)5.0*sC;
  sD=(ValueType)5.0*sD;  

  if (sA>(ValueType)4.0) //1
  {
    result=(sA-(ValueType)4.0)*v5000+sB*v4100+sC*v4010+sD*v4001;
  }
  else if (sB>(ValueType)4.0) //2
  {
    result=sA*v1400+(sB-(ValueType)4.0)*v0500+sC*v0410+sD*v0401;
  }
  else if (sC>(ValueType)4.0) //3
  {
    result=sA*v1040+sB*v0140+(sC-(ValueType)4.0)*v0050+sD*v0041;
  }
  else if (sD>(ValueType)4.0) //4
  {
    result=sA*v1004+sB*v0104+sC*v0014+(sD-(ValueType)4.0)*v0005;
  }

  else if ((sA>(ValueType)3.0)&&(sB>=(ValueType)1.0)) //11
  {
    result=(sA-(ValueType)3.0)*v4100+(sB-(ValueType)1.0)*v3200+sC*v3110+sD*v3101;
  }
  else if ((sB>(ValueType)3.0)&&(sC>=(ValueType)1.0)) //12
  {
    result=sA*v1310+(sB-(ValueType)3.0)*v0410+(sC-(ValueType)1.0)*v0320+sD*v0311;
  }
  else if ((sC>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //13
  {
    result=sA*v1031+sB*v0131+(sC-(ValueType)3.0)*v0041+(sD-(ValueType)1.0)*v0032;
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //14
  {
    result=(sA-(ValueType)1.0)*v2003+sB*v1103+sC*v1013+(sD-(ValueType)3.0)*v1004;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)3.0)) //15
  {
    result=(sA-(ValueType)1.0)*v2300+(sB-(ValueType)3.0)*v1400+sC*v1310+sD*v1301;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)3.0)) //16
  {
    result=sA*v1130+(sB-(ValueType)1.0)*v0230+(sC-(ValueType)3.0)*v0140+sD*v0131;
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //17
  {
    result=sA*v1013+sB*v0113+(sC-(ValueType)1.0)*v0023+(sD-(ValueType)3.0)*v0014;
  }
  else if ((sA>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //18
  {
    result=(sA-(ValueType)3.0)*v4001+sB*v3101+sC*v3011+(sD-(ValueType)1.0)*v3002;
  }
  else if ((sA>(ValueType)3.0)&&(sC>=(ValueType)1.0)) //19
  {
    result=(sA-(ValueType)3.0)*v4010+sB*v3110+(sC-(ValueType)1.0)*v3020+sD*v3011;
  }
  else if ((sB>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //20
  {
    result=sA*v1301+(sB-(ValueType)3.0)*v0401+sC*v0311+(sD-(ValueType)1.0)*v0302;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)3.0)) //21
  {
    result=(sA-(ValueType)1.0)*v2030+sB*v1130+(sC-(ValueType)3.0)*v1040+sD*v1031;
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //22
  {
    result=sA*v1103+(sB-(ValueType)1.0)*v0203+sC*v0113+(sD-(ValueType)3.0)*v0104;
  }

  else if (sA>(ValueType)3.0) //46
  {
    if(axis=='0') //46.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)4.0-sA)*v3110+sD*v4001+s1*v4100+s2*v4010; //46.1.1
        }
        else
        {
          result=sC*v3110+((ValueType)1.0-sB)*v4001+s1*v4100-s2*v3101; //46.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v3110+((ValueType)1.0-sC)*v4001-s1*v3011+s2*v4010; //46.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v3110+(sA-(ValueType)3.0)*v4001-s1*v3011-s2*v3101; //46.1.2
        };
      };
    }
    else if(axis=='1') //46.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v4010+sD*v3101+s1*v4100+s2*v3110; //46.2.1
        }
        else
        {
          result=sC*v4010+((ValueType)4.0-sA)*v3101+s1*v4100-s2*v4001; //46.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)3.0)*v4010+((ValueType)1.0-sC)*v3101-s1*v3011+s2*v3110; //46.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v4010+sB*v3101-s1*v3011-s2*v4001; //46.2.4
        };
      };
    }
    else //46.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v4100+sD*v3011+s1*v3110+s2*v4010; //46.3.1
        }
        else
        {
          result=(sA-(ValueType)3.0)*v4100+((ValueType)1.0-sB)*v3011+s1*v3110-s2*v3101; //46.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v4100+((ValueType)4.0-sA)*v3011-s1*v4001+s2*v4010; //46.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v4100+sC*v3011-s1*v4001-s2*v3101; //46.3.3
        };
      };
    };        
  }
  else if (sB>(ValueType)3.0) //47
  {
    if(axis=='0') //47.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0410+sD*v1301+s1*v1400+s2*v1310; //47.1.1
        }
        else
        {
          result=sC*v0410+((ValueType)4.0-sB)*v1301+s1*v1400-s2*v0401; //47.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)3.0)*v0410+((ValueType)1.0-sC)*v1301-s1*v0311+s2*v1310; //47.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0410+sA*v1301-s1*v0311-s2*v0401; //47.1.2
        };
      };
    }
    else if(axis=='1') //47.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)4.0-sB)*v1310+sD*v0401+s1*v1400+s2*v0410; //47.2.1
        }
        else
        {
          result=sC*v1310+((ValueType)1.0-sA)*v0401+s1*v1400-s2*v1301; //47.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1310+((ValueType)1.0-sC)*v0401-s1*v0311+s2*v0410; //47.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1310+(sB-(ValueType)3.0)*v0401-s1*v0311-s2*v1301; //47.2.4
        };
      };
    }
    else //47.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1400+sD*v0311+s1*v0410+s2*v1310; //47.3.1
        }
        else
        {
          result=sA*v1400+((ValueType)4.0-sB)*v0311+s1*v0410-s2*v0401; //47.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)3.0)*v1400+((ValueType)1.0-sA)*v0311-s1*v1301+s2*v1310; //47.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1400+sC*v0311-s1*v1301-s2*v0401; //47.3.3
        };
      };
    };        
  }
  else if (sC>(ValueType)3.0) //48
  {
    if(axis=='0') //48.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0140+sD*v1031+s1*v1130+s2*v1040; //48.1.1
        }
        else
        {
          result=(sC-(ValueType)3.0)*v0140+((ValueType)1.0-sB)*v1031+s1*v1130-s2*v0131; //48.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0140+((ValueType)4.0-sC)*v1031-s1*v0041+s2*v1040; //48.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0140+sA*v1031-s1*v0041-s2*v0131; //48.1.2
        };
      };
    }
    else if(axis=='1') //48.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1040+sD*v0131+s1*v1130+s2*v0140; //48.2.1
        }
        else
        {
          result=(sC-(ValueType)3.0)*v1040+((ValueType)1.0-sA)*v0131+s1*v1130-s2*v1031; //48.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1040+((ValueType)4.0-sC)*v0131-s1*v0041+s2*v0140; //48.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1040+sB*v0131-s1*v0041-s2*v1031; //48.2.4
        };
      };
    }
    else //48.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)4.0-sC)*v1130+sD*v0041+s1*v0140+s2*v1040; //48.3.1
        }
        else
        {
          result=sA*v1130+((ValueType)1.0-sB)*v0041+s1*v0140-s2*v0131; //48.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1130+((ValueType)1.0-sA)*v0041-s1*v1031+s2*v1040; //48.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1130+(sC-(ValueType)3.0)*v0041-s1*v1031-s2*v0131; //48.3.3
        };
      };
    };        
  }
  else if (sD>(ValueType)3.0) //49
  {
    if(axis=='0') //49.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0113+(sD-(ValueType)3.0)*v1004+s1*v1103+s2*v1013; //49.1.1
        }
        else
        {
          result=sC*v0113+((ValueType)1.0-sB)*v1004+s1*v1103-s2*v0104; //49.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0113+((ValueType)1.0-sC)*v1004-s1*v0014+s2*v1013; //49.1.4
        }
        else
        {
          result=((ValueType)4.0-sD)*v0113+sA*v1004-s1*v0014-s2*v0104; //49.1.2
        };
      };
    }
    else if(axis=='1') //49.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1013+(sD-(ValueType)3.0)*v0104+s1*v1103+s2*v0113; //49.2.1
        }
        else
        {
          result=sC*v1013+((ValueType)1.0-sA)*v0104+s1*v1103-s2*v1004; //49.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1013+((ValueType)1.0-sC)*v0104-s1*v0014+s2*v0113; //49.2.2
        }
        else
        {
          result=((ValueType)4.0-sD)*v1013+sB*v0104-s1*v0014-s2*v1004; //49.2.4
        };
      };
    }
    else //49.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1103+(sD-(ValueType)3.0)*v0014+s1*v0113+s2*v1013; //49.3.1
        }
        else
        {
          result=sA*v1103+((ValueType)1.0-sB)*v0014+s1*v0113-s2*v0104; //49.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1103+((ValueType)1.0-sA)*v0014-s1*v1004+s2*v1013; //49.3.4
        }
        else
        {
          result=((ValueType)4.0-sD)*v1103+sC*v0014-s1*v1004-s2*v0104; //49.3.3
        };
      };
    };        
  }

  else if ((sA>(ValueType)2.0)&&(sB>(ValueType)2.0)) //5
  {
    result=(sA-(ValueType)2.0)*v3200+(sB-(ValueType)2.0)*v2300+sC*v2210+sD*v2201;
  }
  else if ((sB>(ValueType)2.0)&&(sC>(ValueType)2.0)) //6
  {
    result=sA*v1220+(sB-(ValueType)2.0)*v0320+(sC-(ValueType)2.0)*v0230+sD*v0221;
  }
  else if ((sC>(ValueType)2.0)&&(sD>(ValueType)2.0)) //7
  {
    result=sA*v1022+sB*v0122+(sC-(ValueType)2.0)*v0032+(sD-(ValueType)2.0)*v0023;
  }
  else if ((sA>(ValueType)2.0)&&(sD>(ValueType)2.0)) //8
  {
    result=(sA-(ValueType)2.0)*v3002+sB*v2102+sC*v2012+(sD-(ValueType)2.0)*v2003;
  }
  else if ((sA>(ValueType)2.0)&&(sC>(ValueType)2.0)) //9
  {
    result=(sA-(ValueType)2.0)*v3020+sB*v2120+(sC-(ValueType)2.0)*v2030+sD*v2021;
  }
  else if ((sB>(ValueType)2.0)&&(sD>(ValueType)2.0)) //10
  {
    result=sA*v1202+(sB-(ValueType)2.0)*v0302+sC*v0212+(sD-(ValueType)2.0)*v0203;
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)&&(sC>=(ValueType)1.0)) //23
  {
    result=(sA-(ValueType)2.0)*v3110+(sB-(ValueType)1.0)*v2210+(sC-(ValueType)1.0)*v2120+sD*v2111;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //24
  {
    result=(sA-(ValueType)1.0)*v2210+(sB-(ValueType)2.0)*v1310+(sC-(ValueType)1.0)*v1220+sD*v1211;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //25
  {
    result=(sA-(ValueType)1.0)*v2120+(sB-(ValueType)1.0)*v1220+(sC-(ValueType)2.0)*v1130+sD*v1121;
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //26
  {
    result=sA*v1211+(sB-(ValueType)2.0)*v0311+(sC-(ValueType)1.0)*v0221+(sD-(ValueType)1.0)*v0212;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //27
  {
    result=sA*v1121+(sB-(ValueType)1.0)*v0221+(sC-(ValueType)2.0)*v0131+(sD-(ValueType)1.0)*v0122;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //28
  {
    result=sA*v1112+(sB-(ValueType)1.0)*v0212+(sC-(ValueType)1.0)*v0122+(sD-(ValueType)2.0)*v0113;
  }
  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //29
  {
    result=(sA-(ValueType)2.0)*v3101+(sB-(ValueType)1.0)*v2201+sC*v2111+(sD-(ValueType)1.0)*v2102;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //30
  {
    result=(sA-(ValueType)1.0)*v2201+(sB-(ValueType)2.0)*v1301+sC*v1211+(sD-(ValueType)1.0)*v1202;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //31
  {
    result=(sA-(ValueType)1.0)*v2102+(sB-(ValueType)1.0)*v1202+sC*v1112+(sD-(ValueType)2.0)*v1103;
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //32
  {
    result=(sA-(ValueType)2.0)*v3011+sB*v2111+(sC-(ValueType)1.0)*v2021+(sD-(ValueType)1.0)*v2012;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //33
  {
    result=(sA-(ValueType)1.0)*v2021+sB*v1121+(sC-(ValueType)2.0)*v1031+(sD-(ValueType)1.0)*v1022;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //34
  {
    result=(sA-(ValueType)1.0)*v2012+sB*v1112+(sC-(ValueType)1.0)*v1022+(sD-(ValueType)2.0)*v1013;
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)) //50
  {
    if(axis=='0') //50.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sA)*v2210+sD*v3101+s1*v3200+s2*v3110; //50.1.1
        }
        else
        {
          result=sC*v2210+((ValueType)2.0-sB)*v3101+s1*v3200-s2*v2201; //50.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v2210+((ValueType)1.0-sC)*v3101-s1*v2111+s2*v3110; //50.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2210+(sA-(ValueType)2.0)*v3101-s1*v2111-s2*v2201; //50.1.2
        };
      };
    }
    else if(axis=='1') //50.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v3110+sD*v2201+s1*v3200+s2*v2210; //50.2.1
        }
        else
        {
          result=sC*v3110+((ValueType)3.0-sA)*v2201+s1*v3200-s2*v3101; //50.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)2.0)*v3110+((ValueType)1.0-sC)*v2201-s1*v2111+s2*v2210; //50.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v3110+(sB-(ValueType)1.0)*v2201-s1*v2111-s2*v3101; //50.2.4
        };
      };
    }
    else //50.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v3200+sD*v2111+s1*v2210+s2*v3110; //50.3.1
        }
        else
        {
          result=(sA-(ValueType)2.0)*v3200+((ValueType)2.0-sB)*v2111+s1*v2210-s2*v2201; //50.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v3200+((ValueType)3.0-sA)*v2111-s1*v3101+s2*v3110; //50.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v3200+sC*v2111-s1*v3101-s2*v2201; //50.3.3
        };
      };
    };        
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //51
  {
    if(axis=='0') //51.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0320+sD*v1211+s1*v1310+s2*v1220; //51.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v0320+((ValueType)3.0-sB)*v1211+s1*v1310-s2*v0311; //51.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v0320+((ValueType)2.0-sC)*v1211-s1*v0221+s2*v1220; //51.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0320+sA*v1211-s1*v0221-s2*v0311; //51.1.2
        };
      };
    }
    else if(axis=='1') //51.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sB)*v1220+sD*v0311+s1*v1310+s2*v0320; //51.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1220+((ValueType)1.0-sA)*v0311+s1*v1310-s2*v1211; //51.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1220+((ValueType)2.0-sC)*v0311-s1*v0221+s2*v0320; //51.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1220+(sB-(ValueType)2.0)*v0311-s1*v0221-s2*v1211; //51.2.4
        };
      };
    }
    else //51.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v1310+sD*v0221+s1*v0320+s2*v1220; //51.3.1
        }
        else
        {
          result=sA*v1310+((ValueType)3.0-sB)*v0221+s1*v0320-s2*v0311; //51.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v1310+((ValueType)1.0-sA)*v0221-s1*v1211+s2*v1220; //51.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1310+(sC-(ValueType)1.0)*v0221-s1*v1211-s2*v0311; //51.3.3
        };
      };
    };        
  }
  else if ((sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //52
  {
    if(axis=='0') //52.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0131+(sD-(ValueType)1.0)*v1022+s1*v1121+s2*v1031; //52.1.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v0131+((ValueType)1.0-sB)*v1022+s1*v1121-s2*v0122; //52.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0131+((ValueType)3.0-sC)*v1022-s1*v0032+s2*v1031; //52.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v0131+sA*v1022-s1*v0032-s2*v0122; //52.1.2
        };
      };
    }
    else if(axis=='1') //52.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1031+(sD-(ValueType)1.0)*v0122+s1*v1121+s2*v0131; //52.2.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v1031+((ValueType)1.0-sA)*v0122+s1*v1121-s2*v1022; //52.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1031+((ValueType)3.0-sC)*v0122-s1*v0032+s2*v0131; //52.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v1031+sB*v0122-s1*v0032-s2*v1022; //52.2.4
        };
      };
    }
    else //52.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sC)*v1121+(sD-(ValueType)1.0)*v0032+s1*v0131+s2*v1031; //52.3.1
        }
        else
        {
          result=sA*v1121+((ValueType)1.0-sB)*v0032+s1*v0131-s2*v0122; //52.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1121+((ValueType)1.0-sA)*v0032-s1*v1022+s2*v1031; //52.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1121+(sC-(ValueType)2.0)*v0032-s1*v1022-s2*v0122; //52.3.3
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //53
  {
    if(axis=='0') //53.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1112+(sD-(ValueType)2.0)*v2003+s1*v2102+s2*v2012; //53.1.1
        }
        else
        {
          result=sC*v1112+((ValueType)1.0-sB)*v2003+s1*v2102-s2*v1103; //53.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1112+((ValueType)1.0-sC)*v2003-s1*v1013+s2*v2012; //53.1.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v1112+(sA-(ValueType)1.0)*v2003-s1*v1013-s2*v1103; //53.1.2
        };
      };
    }
    else if(axis=='1') //53.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v2012+(sD-(ValueType)2.0)*v1103+s1*v2102+s2*v1112; //53.2.1
        }
        else
        {
          result=sC*v2012+((ValueType)2.0-sA)*v1103+s1*v2102-s2*v2003; //53.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2012+((ValueType)1.0-sC)*v1103-s1*v1013+s2*v1112; //53.2.2
        }
        else
        {
          result=((ValueType)3.0-sD)*v2012+sB*v1103-s1*v1013-s2*v2003; //53.2.4
        };
      };
    }
    else //53.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v2102+(sD-(ValueType)2.0)*v1013+s1*v1112+s2*v2012; //53.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2102+((ValueType)1.0-sB)*v1013+s1*v1112-s2*v1103; //53.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2102+((ValueType)2.0-sA)*v1013-s1*v2003+s2*v2012; //53.3.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v2102+sC*v1013-s1*v2003-s2*v1103; //53.3.3
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)) //54
  {
    if(axis=='0') //54.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1310+sD*v2201+s1*v2300+s2*v2210; //54.1.1
        }
        else
        {
          result=sC*v1310+((ValueType)3.0-sB)*v2201+s1*v2300-s2*v1301; //54.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v1310+((ValueType)1.0-sC)*v2201-s1*v1211+s2*v2210; //54.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1310+(sA-(ValueType)1.0)*v2201-s1*v1211-s2*v1301; //54.1.2
        };
      };
    }
    else if(axis=='1') //54.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sB)*v2210+sD*v1301+s1*v2300+s2*v1310; //54.2.1
        }
        else
        {
          result=sC*v2210+((ValueType)2.0-sA)*v1301+s1*v2300-s2*v2201; //54.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2210+((ValueType)1.0-sC)*v1301-s1*v1211+s2*v1310; //54.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v2210+(sB-(ValueType)2.0)*v1301-s1*v1211-s2*v2201; //54.2.4
        };
      };
    }
    else //54.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v2300+sD*v1211+s1*v1310+s2*v2210; //54.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2300+((ValueType)3.0-sB)*v1211+s1*v1310-s2*v1301; //54.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v2300+((ValueType)2.0-sA)*v1211-s1*v2201+s2*v2210; //54.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2300+sC*v1211-s1*v2201-s2*v1301; //54.3.3
        };
      };
    };        
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //55
  {
    if(axis=='0') //55.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0230+sD*v1121+s1*v1220+s2*v1130; //55.1.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v0230+((ValueType)2.0-sB)*v1121+s1*v1220-s2*v0221; //55.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v0230+((ValueType)3.0-sC)*v1121-s1*v0131+s2*v1130; //55.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v0230+sA*v1121-s1*v0131-s2*v0221; //55.1.2
        };
      };
    }
    else if(axis=='1') //55.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v1130+sD*v0221+s1*v1220+s2*v0230; //55.2.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v1130+((ValueType)1.0-sA)*v0221+s1*v1220-s2*v1121; //55.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1130+((ValueType)3.0-sC)*v0221-s1*v0131+s2*v0230; //55.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v1130+(sB-(ValueType)1.0)*v0221-s1*v0131-s2*v1121; //55.2.4
        };
      };
    }
    else //55.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sC)*v1220+sD*v0131+s1*v0230+s2*v1130; //55.3.1
        }
        else
        {
          result=sA*v1220+((ValueType)2.0-sB)*v0131+s1*v0230-s2*v0221; //55.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1220+((ValueType)1.0-sA)*v0131-s1*v1121+s2*v1130; //55.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1220+(sC-(ValueType)2.0)*v0131-s1*v1121-s2*v0221; //55.3.3
        };
      };
    };        
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //56
  {
    if(axis=='0') //56.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0122+(sD-(ValueType)2.0)*v1013+s1*v1112+s2*v1022; //56.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v0122+((ValueType)1.0-sB)*v1013+s1*v1112-s2*v0113; //56.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v0122+((ValueType)2.0-sC)*v1013-s1*v0023+s2*v1022; //56.1.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v0122+sA*v1013-s1*v0023-s2*v0113; //56.1.2
        };
      };
    }
    else if(axis=='1') //56.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v1022+(sD-(ValueType)2.0)*v0113+s1*v1112+s2*v0122; //56.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1022+((ValueType)1.0-sA)*v0113+s1*v1112-s2*v1013; //56.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1022+((ValueType)2.0-sC)*v0113-s1*v0023+s2*v0122; //56.2.2
        }
        else
        {
          result=((ValueType)3.0-sD)*v1022+sB*v0113-s1*v0023-s2*v1013; //56.2.4
        };
      };
    }
    else //56.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v1112+(sD-(ValueType)2.0)*v0023+s1*v0122+s2*v1022; //56.3.1
        }
        else
        {
          result=sA*v1112+((ValueType)1.0-sB)*v0023+s1*v0122-s2*v0113; //56.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1112+((ValueType)1.0-sA)*v0023-s1*v1013+s2*v1022; //56.3.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v1112+(sC-(ValueType)1.0)*v0023-s1*v1013-s2*v0113; //56.3.3
        };
      };
    };        
  }
  else if ((sA>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //57
  {
    if(axis=='0') //57.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sA)*v2111+(sD-(ValueType)1.0)*v3002+s1*v3101+s2*v3011; //57.1.1
        }
        else
        {
          result=sC*v2111+((ValueType)1.0-sB)*v3002+s1*v3101-s2*v2102; //57.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2111+((ValueType)1.0-sC)*v3002-s1*v2012+s2*v3011; //57.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v2111+(sA-(ValueType)2.0)*v3002-s1*v2012-s2*v2102; //57.1.2
        };
      };
    }
    else if(axis=='1') //57.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v3011+(sD-(ValueType)1.0)*v2102+s1*v3101+s2*v2111; //57.2.1
        }
        else
        {
          result=sC*v3011+((ValueType)3.0-sA)*v2102+s1*v3101-s2*v3002; //57.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)2.0)*v3011+((ValueType)1.0-sC)*v2102-s1*v2012+s2*v2111; //57.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v3011+sB*v2102-s1*v2012-s2*v3002; //57.2.4
        };
      };
    }
    else //57.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v3101+(sD-(ValueType)1.0)*v2012+s1*v2111+s2*v3011; //57.3.1
        }
        else
        {
          result=(sA-(ValueType)2.0)*v3101+((ValueType)1.0-sB)*v2012+s1*v2111-s2*v2102; //57.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v3101+((ValueType)3.0-sA)*v2012-s1*v3002+s2*v3011; //57.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v3101+sC*v2012-s1*v3002-s2*v2102; //57.3.3
        };
      };
    };        
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //58
  {
    if(axis=='0') //58.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sA)*v2120+sD*v3011+s1*v3110+s2*v3020; //58.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v2120+((ValueType)1.0-sB)*v3011+s1*v3110-s2*v2111; //58.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2120+((ValueType)2.0-sC)*v3011-s1*v2021+s2*v3020; //58.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2120+(sA-(ValueType)2.0)*v3011-s1*v2021-s2*v2111; //58.1.2
        };
      };
    }
    else if(axis=='1') //58.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v3020+sD*v2111+s1*v3110+s2*v2120; //58.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v3020+((ValueType)3.0-sA)*v2111+s1*v3110-s2*v3011; //58.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)2.0)*v3020+((ValueType)2.0-sC)*v2111-s1*v2021+s2*v2120; //58.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v3020+sB*v2111-s1*v2021-s2*v3011; //58.2.4
        };
      };
    }
    else //58.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v3110+sD*v2021+s1*v2120+s2*v3020; //58.3.1
        }
        else
        {
          result=(sA-(ValueType)2.0)*v3110+((ValueType)1.0-sB)*v2021+s1*v2120-s2*v2111; //58.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v3110+((ValueType)3.0-sA)*v2021-s1*v3011+s2*v3020; //58.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v3110+(sC-(ValueType)1.0)*v2021-s1*v3011-s2*v2111; //58.3.3
        };
      };
    };        
  }
  else if ((sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //59
  {
    if(axis=='0') //59.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0311+(sD-(ValueType)1.0)*v1202+s1*v1301+s2*v1211; //59.1.1
        }
        else
        {
          result=sC*v0311+((ValueType)3.0-sB)*v1202+s1*v1301-s2*v0302; //59.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v0311+((ValueType)1.0-sC)*v1202-s1*v0212+s2*v1211; //59.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v0311+sA*v1202-s1*v0212-s2*v0302; //59.1.2
        };
      };
    }
    else if(axis=='1') //59.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sB)*v1211+(sD-(ValueType)1.0)*v0302+s1*v1301+s2*v0311; //59.2.1
        }
        else
        {
          result=sC*v1211+((ValueType)1.0-sA)*v0302+s1*v1301-s2*v1202; //59.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1211+((ValueType)1.0-sC)*v0302-s1*v0212+s2*v0311; //59.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v1211+(sB-(ValueType)2.0)*v0302-s1*v0212-s2*v1202; //59.2.4
        };
      };
    }
    else //59.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1301+(sD-(ValueType)1.0)*v0212+s1*v0311+s2*v1211; //59.3.1
        }
        else
        {
          result=sA*v1301+((ValueType)3.0-sB)*v0212+s1*v0311-s2*v0302; //59.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)2.0)*v1301+((ValueType)1.0-sA)*v0212-s1*v1202+s2*v1211; //59.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1301+sC*v0212-s1*v1202-s2*v0302; //59.3.3
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //60
  {
    if(axis=='0') //60.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1130+sD*v2021+s1*v2120+s2*v2030; //60.1.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v1130+((ValueType)1.0-sB)*v2021+s1*v2120-s2*v1121; //60.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1130+((ValueType)3.0-sC)*v2021-s1*v1031+s2*v2030; //60.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1130+(sA-(ValueType)1.0)*v2021-s1*v1031-s2*v1121; //60.1.2
        };
      };
    }
    else if(axis=='1') //60.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v2030+sD*v1121+s1*v2120+s2*v1130; //60.2.1
        }
        else
        {
          result=(sC-(ValueType)2.0)*v2030+((ValueType)2.0-sA)*v1121+s1*v2120-s2*v2021; //60.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2030+((ValueType)3.0-sC)*v1121-s1*v1031+s2*v1130; //60.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v2030+sB*v1121-s1*v1031-s2*v2021; //60.2.4
        };
      };
    }
    else //60.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)3.0-sC)*v2120+sD*v1031+s1*v1130+s2*v2030; //60.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2120+((ValueType)1.0-sB)*v1031+s1*v1130-s2*v1121; //60.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2120+((ValueType)2.0-sA)*v1031-s1*v2021+s2*v2030; //60.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2120+(sC-(ValueType)2.0)*v1031-s1*v2021-s2*v1121; //60.3.3
        };
      };
    };        
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //61
  {
    if(axis=='0') //61.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0212+(sD-(ValueType)2.0)*v1103+s1*v1202+s2*v1112; //61.1.1
        }
        else
        {
          result=sC*v0212+((ValueType)2.0-sB)*v1103+s1*v1202-s2*v0203; //61.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v0212+((ValueType)1.0-sC)*v1103-s1*v0113+s2*v1112; //61.1.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v0212+sA*v1103-s1*v0113-s2*v0203; //61.1.2
        };
      };
    }
    else if(axis=='1') //61.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v1112+(sD-(ValueType)2.0)*v0203+s1*v1202+s2*v0212; //61.2.1
        }
        else
        {
          result=sC*v1112+((ValueType)1.0-sA)*v0203+s1*v1202-s2*v1103; //61.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1112+((ValueType)1.0-sC)*v0203-s1*v0113+s2*v0212; //61.2.2
        }
        else
        {
          result=((ValueType)3.0-sD)*v1112+(sB-(ValueType)1.0)*v0203-s1*v0113-s2*v1103; //61.2.4
        };
      };
    }
    else //61.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v1202+(sD-(ValueType)2.0)*v0113+s1*v0212+s2*v1112; //61.3.1
        }
        else
        {
          result=sA*v1202+((ValueType)2.0-sB)*v0113+s1*v0212-s2*v0203; //61.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1202+((ValueType)1.0-sA)*v0113-s1*v1103+s2*v1112; //61.3.4
        }
        else
        {
          result=((ValueType)3.0-sD)*v1202+sC*v0113-s1*v1103-s2*v0203; //61.3.3
        };
      };
    };        
  }

  else if (sA>(ValueType)2.0) //36
  {
    result=((ValueType)3.0-sA)*v2111+((ValueType)1.0-sB)*v3011+((ValueType)1.0-sC)*v3101+((ValueType)1.0-sD)*v3110;
  }
  else if (sB>(ValueType)2.0) //37
  {
    result=((ValueType)1.0-sA)*v0311+((ValueType)3.0-sB)*v1211+((ValueType)1.0-sC)*v1301+((ValueType)1.0-sD)*v1310;
  }
  else if (sC>(ValueType)2.0) //38
  {
    result=((ValueType)1.0-sA)*v0131+((ValueType)1.0-sB)*v1031+((ValueType)3.0-sC)*v1121+((ValueType)1.0-sD)*v1130;
  }
  else if (sD>(ValueType)2.0) //39
  {
    result=((ValueType)1.0-sA)*v0113+((ValueType)1.0-sB)*v1013+((ValueType)1.0-sC)*v1103+((ValueType)3.0-sD)*v1112;
  }

  else if ((sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //40
  {
    result=((ValueType)2.0-sA)*v1211+((ValueType)2.0-sB)*v2111+((ValueType)1.0-sC)*v2201+((ValueType)1.0-sD)*v2210;
  }
  else if ((sA<(ValueType)1.0)&&(sD<(ValueType)1.0)) //41
  {
    result=((ValueType)1.0-sA)*v0221+((ValueType)2.0-sB)*v1121+((ValueType)2.0-sC)*v1211+((ValueType)1.0-sD)*v1220;
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //42
  {
    result=((ValueType)1.0-sA)*v0122+((ValueType)1.0-sB)*v1022+((ValueType)2.0-sC)*v1112+((ValueType)2.0-sD)*v1121;
  }
  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //43
  {
    result=((ValueType)2.0-sA)*v1112+((ValueType)1.0-sB)*v2012+((ValueType)1.0-sC)*v2102+((ValueType)2.0-sD)*v2111;
  }
  else if ((sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //44
  {
    result=((ValueType)2.0-sA)*v1121+((ValueType)1.0-sB)*v2021+((ValueType)2.0-sC)*v2111+((ValueType)1.0-sD)*v2120;
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //45
  {
    result=((ValueType)1.0-sA)*v0212+((ValueType)2.0-sB)*v1112+((ValueType)1.0-sC)*v1202+((ValueType)2.0-sD)*v1211;
  }

  else if (sA<(ValueType)1.0) //62
  {
    if(axis=='0') //62.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sA)*v0221+(sD-(ValueType)1.0)*v1112+s1*v1211+s2*v1121; //62.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v0221+((ValueType)2.0-sB)*v1112+s1*v1211-s2*v0212; //62.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v0221+((ValueType)2.0-sC)*v1112-s1*v0122+s2*v1121; //62.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v0221+sA*v1112-s1*v0122-s2*v0212; //62.1.2
        };
      };
    }
    else if(axis=='1') //62.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v1121+(sD-(ValueType)1.0)*v0212+s1*v1211+s2*v0221; //62.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1121+((ValueType)1.0-sA)*v0212+s1*v1211-s2*v1112; //62.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sA*v1121+((ValueType)2.0-sC)*v0212-s1*v0122+s2*v0221; //62.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v1121+(sB-(ValueType)1.0)*v0212-s1*v0122-s2*v1112; //62.2.4
        };
      };
    }
    else //62.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v1211+(sD-(ValueType)1.0)*v0122+s1*v0221+s2*v1121; //62.3.1
        }
        else
        {
          result=sA*v1211+((ValueType)2.0-sB)*v0122+s1*v0221-s2*v0212; //62.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1211+((ValueType)1.0-sA)*v0122-s1*v1112+s2*v1121; //62.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1211+(sC-(ValueType)1.0)*v0122-s1*v1112-s2*v0212; //62.3.3
        };
      };
    };        
  }
  else if (sB<(ValueType)1.0) //63
  {
    if(axis=='0') //63.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1121+(sD-(ValueType)1.0)*v2012+s1*v2111+s2*v2021; //63.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1121+((ValueType)1.0-sB)*v2012+s1*v2111-s2*v1112; //63.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v1121+((ValueType)2.0-sC)*v2012-s1*v1022+s2*v2021; //63.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1121+(sA-(ValueType)1.0)*v2012-s1*v1022-s2*v1112; //63.1.2
        };
      };
    }
    else if(axis=='1') //63.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sB)*v2021+(sD-(ValueType)1.0)*v1112+s1*v2111+s2*v1121; //63.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v2021+((ValueType)2.0-sA)*v1112+s1*v2111-s2*v2012; //63.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2021+((ValueType)2.0-sC)*v1112-s1*v1022+s2*v1121; //63.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v2021+sB*v1112-s1*v1022-s2*v2012; //63.2.4
        };
      };
    }
    else //63.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v2111+(sD-(ValueType)1.0)*v1022+s1*v1121+s2*v2021; //63.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2111+((ValueType)1.0-sB)*v1022+s1*v1121-s2*v1112; //63.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=sB*v2111+((ValueType)2.0-sA)*v1022-s1*v2012+s2*v2021; //63.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v2111+(sC-(ValueType)1.0)*v1022-s1*v2012-s2*v1112; //63.3.3
        };
      };
    };        
  }
  else if (sC<(ValueType)1.0) //64
  {
    if(axis=='0') //64.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1211+(sD-(ValueType)1.0)*v2102+s1*v2201+s2*v2111; //64.1.1
        }
        else
        {
          result=sC*v1211+((ValueType)2.0-sB)*v2102+s1*v2201-s2*v1202; //64.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1211+((ValueType)1.0-sC)*v2102-s1*v1112+s2*v2111; //64.1.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v1211+(sA-(ValueType)1.0)*v2102-s1*v1112-s2*v1202; //64.1.2
        };
      };
    }
    else if(axis=='1') //64.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v2111+(sD-(ValueType)1.0)*v1202+s1*v2201+s2*v1211; //64.2.1
        }
        else
        {
          result=sC*v2111+((ValueType)2.0-sA)*v1202+s1*v2201-s2*v2102; //64.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2111+((ValueType)1.0-sC)*v1202-s1*v1112+s2*v1211; //64.2.2
        }
        else
        {
          result=((ValueType)2.0-sD)*v2111+(sB-(ValueType)1.0)*v1202-s1*v1112-s2*v2102; //64.2.4
        };
      };
    }
    else //64.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)1.0-sC)*v2201+(sD-(ValueType)1.0)*v1112+s1*v1211+s2*v2111; //64.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2201+((ValueType)2.0-sB)*v1112+s1*v1211-s2*v1202; //64.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v2201+((ValueType)2.0-sA)*v1112-s1*v2102+s2*v2111; //64.3.4
        }
        else
        {
          result=((ValueType)2.0-sD)*v2201+sC*v1112-s1*v2102-s2*v1202; //64.3.3
        };
      };
    };        
  }
  else if (sD<(ValueType)1.0) //65
  {
    if(axis=='0') //65.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sA)*v1220+sD*v2111+s1*v2210+s2*v2120; //65.1.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v1220+((ValueType)2.0-sB)*v2111+s1*v2210-s2*v1211; //65.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v1220+((ValueType)2.0-sC)*v2111-s1*v1121+s2*v2120; //65.1.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v1220+(sA-(ValueType)1.0)*v2111-s1*v1121-s2*v1211; //65.1.2
        };
      };
    }
    else if(axis=='1') //65.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sB)*v2120+sD*v1211+s1*v2210+s2*v1220; //65.2.1
        }
        else
        {
          result=(sC-(ValueType)1.0)*v2120+((ValueType)2.0-sA)*v1211+s1*v2210-s2*v2111; //65.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sA-(ValueType)1.0)*v2120+((ValueType)2.0-sC)*v1211-s1*v1121+s2*v1220; //65.2.2
        }
        else
        {
          result=((ValueType)1.0-sD)*v2120+(sB-(ValueType)1.0)*v1211-s1*v1121-s2*v2111; //65.2.4
        };
      };
    }
    else //65.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          result=((ValueType)2.0-sC)*v2210+sD*v1121+s1*v1220+s2*v2120; //65.3.1
        }
        else
        {
          result=(sA-(ValueType)1.0)*v2210+((ValueType)2.0-sB)*v1121+s1*v1220-s2*v1211; //65.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          result=(sB-(ValueType)1.0)*v2210+((ValueType)2.0-sA)*v1121-s1*v2111+s2*v2120; //65.3.4
        }
        else
        {
          result=((ValueType)1.0-sD)*v2210+(sC-(ValueType)1.0)*v1121-s1*v2111-s2*v1211; //65.3.3
        };
      };
    };        
  }

  else //35
  {
    result=(sA-(ValueType)1.0)*v2111+(sB-(ValueType)1.0)*v1211+(sC-(ValueType)1.0)*v1121+(sD-(ValueType)1.0)*v1112;
  };

  return(result);
}

template <typename ValueType>
ValueType Interpolate_3D_Mono_5(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v5000, const ValueType v0500, const ValueType v0050, const ValueType v0005,
const ValueType v4100, const ValueType v0410, const ValueType v0041, const ValueType v1004, 
const ValueType v1400, const ValueType v0140, const ValueType v0014, const ValueType v4001, 
const ValueType v4010, const ValueType v0401, const ValueType v1040, const ValueType v0104,
const ValueType v3200, const ValueType v0320, const ValueType v0032, const ValueType v2003, 
const ValueType v2300, const ValueType v0230, const ValueType v0023, const ValueType v3002, 
const ValueType v3020, const ValueType v0302, const ValueType v2030, const ValueType v0203,
const ValueType v3110, const ValueType v0311, const ValueType v1031, const ValueType v1103, 
const ValueType v1310, const ValueType v0131, const ValueType v1013, const ValueType v3101, 
const ValueType v1130, const ValueType v0113, const ValueType v3011, const ValueType v1301,
const ValueType v2210, const ValueType v0221, const ValueType v1022, const ValueType v2102,
const ValueType v1220, const ValueType v0122, const ValueType v2012, const ValueType v2201, 
const ValueType v2120, const ValueType v0212, const ValueType v2021, const ValueType v1202, 
const ValueType v2111, const ValueType v1211, const ValueType v1121, const ValueType v1112)
{
  ValueType result;
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;
  ValueType Vmin, Vmax;
  ValueType w5000, w0500, w0050, w0005;
  ValueType w4100, w0410, w0041, w1004; 
  ValueType w1400, w0140, w0014, w4001; 
  ValueType w4010, w0401, w1040, w0104;
  ValueType w3200, w0320, w0032, w2003; 
  ValueType w2300, w0230, w0023, w3002; 
  ValueType w3020, w0302, w2030, w0203;
  ValueType w3110, w0311, w1031, w1103; 
  ValueType w1310, w0131, w1013, w3101; 
  ValueType w1130, w0113, w3011, w1301;
  ValueType w2210, w0221, w1022, w2102;
  ValueType w1220, w0122, w2012, w2201; 
  ValueType w2120, w0212, w2021, w1202; 
  ValueType w2111, w1211, w1121, w1112;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  w5000=sA*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*((ValueType)5.0*sA-(ValueType)4.0)/(ValueType)24.0;
  w0500=sB*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*((ValueType)5.0*sB-(ValueType)4.0)/(ValueType)24.0;
  w0050=sC*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*((ValueType)5.0*sC-(ValueType)4.0)/(ValueType)24.0;
  w0005=sD*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*((ValueType)5.0*sD-(ValueType)4.0)/(ValueType)24.0;

  w4100=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0410=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0041=sC*sD*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w1004=sD*sA*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w1400=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0140=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0014=sD*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w4001=sA*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w4010=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*((ValueType)5.0*sA-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0401=sB*sD*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*((ValueType)5.0*sB-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w1040=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*((ValueType)5.0*sC-(ValueType)3.0)*(ValueType)(25.0/24.0);
  w0104=sD*sB*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*((ValueType)5.0*sD-(ValueType)3.0)*(ValueType)(25.0/24.0);

  w3200=sA*sB*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0320=sB*sC*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0032=sC*sD*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w2003=sD*sA*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w2300=sB*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0230=sC*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0023=sD*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w3002=sA*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w3020=sA*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0302=sB*sD*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w2030=sC*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(25.0/12.0);
  w0203=sD*sB*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(25.0/12.0);

  w3110=sA*sB*sC*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w0311=sB*sC*sD*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1031=sC*sD*sA*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1103=sD*sA*sB*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1310=sB*sC*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w0131=sC*sD*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1013=sD*sA*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w3101=sA*sB*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1130=sC*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w0113=sD*sB*sC*((ValueType)5.0*sD-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w3011=sA*sC*sD*((ValueType)5.0*sA-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)2.0)*(ValueType)(125.0/6.0);
  w1301=sB*sD*sA*((ValueType)5.0*sB-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)2.0)*(ValueType)(125.0/6.0);

  w2210=sA*sB*sC*(5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w0221=sB*sC*sD*(5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w1022=sC*sD*sA*(5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2102=sD*sA*sB*(5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w1220=sB*sC*sA*(5.0*sB-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w0122=sC*sD*sB*(5.0*sC-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2012=sD*sA*sC*(5.0*sD-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2201=sA*sB*sD*(5.0*sA-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2120=sC*sA*sB*(5.0*sC-(ValueType)1.0)*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w0212=sD*sB*sC*(5.0*sD-(ValueType)1.0)*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w2021=sA*sC*sD*(5.0*sA-(ValueType)1.0)*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(125.0/4.0);
  w1202=sB*sD*sA*(5.0*sB-(ValueType)1.0)*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(125.0/4.0);

  w2111=sA*sB*sC*sD*((ValueType)5.0*sA-(ValueType)1.0)*(ValueType)(625.0/2.0);
  w1211=sB*sC*sD*sA*((ValueType)5.0*sB-(ValueType)1.0)*(ValueType)(625.0/2.0);
  w1121=sC*sD*sA*sB*((ValueType)5.0*sC-(ValueType)1.0)*(ValueType)(625.0/2.0);
  w1112=sD*sA*sB*sC*((ValueType)5.0*sD-(ValueType)1.0)*(ValueType)(625.0/2.0);

  result=v5000*w5000+v0500*w0500+v0050*w0050+v0005*w0005+
         v4100*w4100+v0410*w0410+v0041*w0041+v1004*w1004+
         v1400*w1400+v0140*w0140+v0014*w0014+v4001*w4001+
         v4010*w4010+v0401*w0401+v1040*w1040+v0104*w0104+
         v3200*w3200+v0320*w0320+v0032*w0032+v2003*w2003+
         v2300*w2300+v0230*w0230+v0023*w0023+v3002*w3002+
         v3020*w3020+v0302*w0302+v2030*w2030+v0203*w0203+
         v3110*w3110+v0311*w0311+v1031*w1031+v1103*w1103+
         v1310*w1310+v0131*w0131+v1013*w1013+v3101*w3101+
         v1130*w1130+v0113*w0113+v3011*w3011+v1301*w1301+
         v2210*w2210+v0221*w0221+v1022*w1022+v2102*w2102+
         v1220*w1220+v0122*w0122+v2012*w2012+v2201*w2201+
         v2120*w2120+v0212*w0212+v2021*w2021+v1202*w1202+
         v2111*w2111+v1211*w1211+v1121*w1121+v1112*w1112;

  sA=(ValueType)5.0*sA;
  sB=(ValueType)5.0*sB;
  sC=(ValueType)5.0*sC;
  sD=(ValueType)5.0*sD;  

  if (sA>(ValueType)4.0) //1
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v5000, v4100, v4010, v4001);
  }
  else if (sB>(ValueType)4.0) //2
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1400, v0500, v0410, v0401);
  }
  else if (sC>(ValueType)4.0) //3
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1040, v0140, v0050, v0041);
  }
  else if (sD>(ValueType)4.0) //4
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1004, v0104, v0014, v0005);
  }

  else if ((sA>(ValueType)3.0)&&(sB>=(ValueType)1.0)) //11
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4100, v3200, v3110, v3101);
  }
  else if ((sB>(ValueType)3.0)&&(sC>=(ValueType)1.0)) //12
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0410, v0320, v0311);
  }
  else if ((sC>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //13
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1031, v0131, v0041, v0032);
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //14
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2003, v1103, v1013, v1004);
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)3.0)) //15
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2300, v1400, v1310, v1301);
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)3.0)) //16
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0230, v0140, v0131);
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //17
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1013, v0113, v0023, v0014);
  }
  else if ((sA>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //18
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4001, v3101, v3011, v3002);
  }
  else if ((sA>(ValueType)3.0)&&(sC>=(ValueType)1.0)) //19
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4010, v3110, v3020, v3011);
  }
  else if ((sB>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //20
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1301, v0401, v0311, v0302);
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)3.0)) //21
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2030, v1130, v1040, v1031);
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //22
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1103, v0203, v0113, v0104);
  }

  else if (sA>(ValueType)3.0) //46
  {
    if(axis=='0') //46.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v4001, v4100, v4010); //46.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v4001, v4100, v3101); //46.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v4001, v3011, v4010); //46.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v4001, v3011, v3101); //46.1.2
        };
      };
    }
    else if(axis=='1') //46.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4010, v3101, v4100, v3110); //46.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4010, v3101, v4100, v4001); //46.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4010, v3101, v3011, v3110); //46.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4010, v3101, v3011, v4001); //46.2.4
        };
      };
    }
    else //46.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4100, v3011, v3110, v4010); //46.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4100, v3011, v3110, v3101); //46.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4100, v3011, v4001, v4010); //46.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v4100, v3011, v4001, v3101); //46.3.3
        };
      };
    };        
  }
  else if (sB>(ValueType)3.0) //47
  {
    if(axis=='0') //47.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0410, v1301, v1400, v1310); //47.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0410, v1301, v1400, v0401); //47.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0410, v1301, v0311, v1310); //47.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0410, v1301, v0311, v0401); //47.1.2
        };
      };
    }
    else if(axis=='1') //47.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0401, v1400, v0410); //47.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0401, v1400, v1301); //47.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0401, v0311, v0410); //47.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0401, v0311, v1301); //47.2.4
        };
      };
    }
    else //47.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1400, v0311, v0410, v1310); //47.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1400, v0311, v0410, v0401); //47.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1400, v0311, v1301, v1310); //47.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1400, v0311, v1301, v0401); //47.3.3
        };
      };
    };        
  }
  else if (sC>(ValueType)3.0) //48
  {
    if(axis=='0') //48.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0140, v1031, v1130, v1040); //48.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0140, v1031, v1130, v0131); //48.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0140, v1031, v0041, v1040); //48.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0140, v1031, v0041, v0131); //48.1.2
        };
      };
    }
    else if(axis=='1') //48.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1040, v0131, v1130, v0140); //48.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1040, v0131, v1130, v1031); //48.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1040, v0131, v0041, v0140); //48.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1040, v0131, v0041, v1031); //48.2.4
        };
      };
    }
    else //48.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0041, v0140, v1040); //48.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0041, v0140, v0131); //48.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0041, v1031, v1040); //48.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0041, v1031, v0131); //48.3.3
        };
      };
    };        
  }
  else if (sD>(ValueType)3.0) //49
  {
    if(axis=='0') //49.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0113, v1004, v1103, v1013); //49.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0113, v1004, v1103, v0104); //49.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0113, v1004, v0014, v1013); //49.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0113, v1004, v0014, v0104); //49.1.2
        };
      };
    }
    else if(axis=='1') //49.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1013, v0104, v1103, v0113); //49.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1013, v0104, v1103, v1004); //49.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1013, v0104, v0014, v0113); //49.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1013, v0104, v0014, v1004); //49.2.4
        };
      };
    }
    else //49.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1103, v0014, v0113, v1013); //49.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1103, v0014, v0113, v0104); //49.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1103, v0014, v1004, v1013); //49.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1103, v0014, v1004, v0104); //49.3.3
        };
      };
    };        
  }

  else if ((sA>(ValueType)2.0)&&(sB>(ValueType)2.0)) //5
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3200, v2300, v2210, v2201);
  }
  else if ((sB>(ValueType)2.0)&&(sC>(ValueType)2.0)) //6
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0320, v0230, v0221);
  }
  else if ((sC>(ValueType)2.0)&&(sD>(ValueType)2.0)) //7
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1022, v0122, v0032, v0023);
  }
  else if ((sA>(ValueType)2.0)&&(sD>(ValueType)2.0)) //8
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3002, v2102, v2012, v2003);
  }
  else if ((sA>(ValueType)2.0)&&(sC>(ValueType)2.0)) //9
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3020, v2120, v2030, v2021);
  }
  else if ((sB>(ValueType)2.0)&&(sD>(ValueType)2.0)) //10
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1202, v0302, v0212, v0203);
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)&&(sC>=(ValueType)1.0)) //23
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2210, v2120, v2111);
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //24
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1310, v1220, v1211);
  }
  else if ((sA>=(ValueType)1.0)&&(sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //25
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1220, v1130, v1121);
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //26
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0311, v0221, v0212);
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //27
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0221, v0131, v0122);
  }
  else if ((sB>=(ValueType)1.0)&&(sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //28
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0212, v0122, v0113);
  }
  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //29
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3101, v2201, v2111, v2102);
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //30
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2201, v1301, v1211, v1202);
  }
  else if ((sA>=(ValueType)1.0)&&(sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //31
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2102, v1202, v1112, v1103);
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //32
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3011, v2111, v2021, v2012);
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //33
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2021, v1121, v1031, v1022);
  }
  else if ((sA>=(ValueType)1.0)&&(sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //34
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2012, v1112, v1022, v1013);
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)) //50
  {
    if(axis=='0') //50.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v3101, v3200, v3110); //50.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v3101, v3200, v2201); //50.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v3101, v2111, v3110); //50.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v3101, v2111, v2201); //50.1.2
        };
      };
    }
    else if(axis=='1') //50.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2201, v3200, v2210); //50.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2201, v3200, v3101); //50.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2201, v2111, v2210); //50.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2201, v2111, v3101); //50.2.4
        };
      };
    }
    else //50.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3200, v2111, v2210, v3110); //50.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3200, v2111, v2210, v2201); //50.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3200, v2111, v3101, v3110); //50.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3200, v2111, v3101, v2201); //50.3.3
        };
      };
    };        
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //51
  {
    if(axis=='0') //51.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0320, v1211, v1310, v1220); //51.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0320, v1211, v1310, v0311); //51.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0320, v1211, v0221, v1220); //51.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0320, v1211, v0221, v0311); //51.1.2
        };
      };
    }
    else if(axis=='1') //51.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0311, v1310, v0320); //51.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0311, v1310, v1211); //51.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0311, v0221, v0320); //51.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0311, v0221, v1211); //51.2.4
        };
      };
    }
    else //51.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0221, v0320, v1220); //51.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0221, v0320, v0311); //51.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0221, v1211, v1220); //51.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v0221, v1211, v0311); //51.3.3
        };
      };
    };        
  }
  else if ((sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //52
  {
    if(axis=='0') //52.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0131, v1022, v1121, v1031); //52.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0131, v1022, v1121, v0122); //52.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0131, v1022, v0032, v1031); //52.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0131, v1022, v0032, v0122); //52.1.2
        };
      };
    }
    else if(axis=='1') //52.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1031, v0122, v1121, v0131); //52.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1031, v0122, v1121, v1022); //52.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1031, v0122, v0032, v0131); //52.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1031, v0122, v0032, v1022); //52.2.4
        };
      };
    }
    else //52.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0032, v0131, v1031); //52.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0032, v0131, v0122); //52.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0032, v1022, v1031); //52.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0032, v1022, v0122); //52.3.3
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //53
  {
    if(axis=='0') //53.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v2003, v2102, v2012); //53.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v2003, v2102, v1103); //53.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v2003, v1013, v2012); //53.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v2003, v1013, v1103); //53.1.2
        };
      };
    }
    else if(axis=='1') //53.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2012, v1103, v2102, v1112); //53.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2012, v1103, v2102, v2003); //53.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2012, v1103, v1013, v1112); //53.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2012, v1103, v1013, v2003); //53.2.4
        };
      };
    }
    else //53.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2102, v1013, v1112, v2012); //53.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2102, v1013, v1112, v1103); //53.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2102, v1013, v2003, v2012); //53.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2102, v1013, v2003, v1103); //53.3.3
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)) //54
  {
    if(axis=='0') //54.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v2201, v2300, v2210); //54.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v2201, v2300, v1301); //54.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v2201, v1211, v2210); //54.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1310, v2201, v1211, v1301); //54.1.2
        };
      };
    }
    else if(axis=='1') //54.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1301, v2300, v1310); //54.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1301, v2300, v2201); //54.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1301, v1211, v1310); //54.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1301, v1211, v2201); //54.2.4
        };
      };
    }
    else //54.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2300, v1211, v1310, v2210); //54.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2300, v1211, v1310, v1301); //54.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2300, v1211, v2201, v2210); //54.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2300, v1211, v2201, v1301); //54.3.3
        };
      };
    };        
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //55
  {
    if(axis=='0') //55.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0230, v1121, v1220, v1130); //55.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0230, v1121, v1220, v0221); //55.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0230, v1121, v0131, v1130); //55.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0230, v1121, v0131, v0221); //55.1.2
        };
      };
    }
    else if(axis=='1') //55.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0221, v1220, v0230); //55.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0221, v1220, v1121); //55.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0221, v0131, v0230); //55.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v0221, v0131, v1121); //55.2.4
        };
      };
    }
    else //55.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0131, v0230, v1130); //55.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0131, v0230, v0221); //55.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0131, v1121, v1130); //55.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v0131, v1121, v0221); //55.3.3
        };
      };
    };        
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //56
  {
    if(axis=='0') //56.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0122, v1013, v1112, v1022); //56.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0122, v1013, v1112, v0113); //56.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0122, v1013, v0023, v1022); //56.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0122, v1013, v0023, v0113); //56.1.2
        };
      };
    }
    else if(axis=='1') //56.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1022, v0113, v1112, v0122); //56.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1022, v0113, v1112, v1013); //56.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1022, v0113, v0023, v0122); //56.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1022, v0113, v0023, v1013); //56.2.4
        };
      };
    }
    else //56.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0023, v0122, v1022); //56.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0023, v0122, v0113); //56.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0023, v1013, v1022); //56.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0023, v1013, v0113); //56.3.3
        };
      };
    };        
  }
  else if ((sA>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //57
  {
    if(axis=='0') //57.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v3002, v3101, v3011); //57.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v3002, v3101, v2102); //57.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v3002, v2012, v3011); //57.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v3002, v2012, v2102); //57.1.2
        };
      };
    }
    else if(axis=='1') //57.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3011, v2102, v3101, v2111); //57.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3011, v2102, v3101, v3002); //57.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3011, v2102, v2012, v2111); //57.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3011, v2102, v2012, v3002); //57.2.4
        };
      };
    }
    else //57.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3101, v2012, v2111, v3011); //57.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3101, v2012, v2111, v2102); //57.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3101, v2012, v3002, v3011); //57.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3101, v2012, v3002, v2102); //57.3.3
        };
      };
    };        
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //58
  {
    if(axis=='0') //58.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v3011, v3110, v3020); //58.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v3011, v3110, v2111); //58.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v3011, v2021, v3020); //58.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v3011, v2021, v2111); //58.1.2
        };
      };
    }
    else if(axis=='1') //58.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3020, v2111, v3110, v2120); //58.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3020, v2111, v3110, v3011); //58.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3020, v2111, v2021, v2120); //58.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3020, v2111, v2021, v3011); //58.2.4
        };
      };
    }
    else //58.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2021, v2120, v3020); //58.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2021, v2120, v2111); //58.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2021, v3011, v3020); //58.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v3110, v2021, v3011, v2111); //58.3.3
        };
      };
    };        
  }
  else if ((sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //59
  {
    if(axis=='0') //59.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0311, v1202, v1301, v1211); //59.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0311, v1202, v1301, v0302); //59.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0311, v1202, v0212, v1211); //59.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0311, v1202, v0212, v0302); //59.1.2
        };
      };
    }
    else if(axis=='1') //59.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0302, v1301, v0311); //59.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0302, v1301, v1202); //59.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0302, v0212, v0311); //59.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0302, v0212, v1202); //59.2.4
        };
      };
    }
    else //59.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1301, v0212, v0311, v1211); //59.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1301, v0212, v0311, v0302); //59.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1301, v0212, v1202, v1211); //59.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1301, v0212, v1202, v0302); //59.3.3
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //60
  {
    if(axis=='0') //60.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v2021, v2120, v2030); //60.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v2021, v2120, v1121); //60.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v2021, v1031, v2030); //60.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1130, v2021, v1031, v1121); //60.1.2
        };
      };
    }
    else if(axis=='1') //60.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2030, v1121, v2120, v1130); //60.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2030, v1121, v2120, v2021); //60.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2030, v1121, v1031, v1130); //60.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2030, v1121, v1031, v2021); //60.2.4
        };
      };
    }
    else //60.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1031, v1130, v2030); //60.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1031, v1130, v1121); //60.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1031, v2021, v2030); //60.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1031, v2021, v1121); //60.3.3
        };
      };
    };        
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //61
  {
    if(axis=='0') //61.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0212, v1103, v1202, v1112); //61.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0212, v1103, v1202, v0203); //61.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0212, v1103, v0113, v1112); //61.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0212, v1103, v0113, v0203); //61.1.2
        };
      };
    }
    else if(axis=='1') //61.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0203, v1202, v0212); //61.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0203, v1202, v1103); //61.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0203, v0113, v0212); //61.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v0203, v0113, v1103); //61.2.4
        };
      };
    }
    else //61.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1202, v0113, v0212, v1112); //61.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1202, v0113, v0212, v0203); //61.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1202, v0113, v1103, v1112); //61.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1202, v0113, v1103, v0203); //61.3.3
        };
      };
    };        
  }

  else if (sA>(ValueType)2.0) //36
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v3011, v3101, v3110);
  }
  else if (sB>(ValueType)2.0) //37
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0311, v1211, v1301, v1310);
  }
  else if (sC>(ValueType)2.0) //38
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0131, v1031, v1121, v1130);
  }
  else if (sD>(ValueType)2.0) //39
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0113, v1013, v1103, v1112);
  }

  else if ((sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //40
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v2111, v2201, v2210);
  }
  else if ((sA<(ValueType)1.0)&&(sD<(ValueType)1.0)) //41
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0221, v1121, v1211, v1220);
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //42
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0122, v1022, v1112, v1121);
  }
  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //43
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1112, v2012, v2102, v2111);
  }
  else if ((sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //44
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v2021, v2111, v2120);
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //45
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0212, v1112, v1202, v1211);
  }

  else if (sA<(ValueType)1.0) //62
  {
    if(axis=='0') //62.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0221, v1112, v1211, v1121); //62.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0221, v1112, v1211, v0212); //62.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0221, v1112, v0122, v1121); //62.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v0221, v1112, v0122, v0212); //62.1.2
        };
      };
    }
    else if(axis=='1') //62.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0212, v1211, v0221); //62.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0212, v1211, v1112); //62.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0212, v0122, v0221); //62.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v0212, v0122, v1112); //62.2.4
        };
      };
    }
    else //62.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0122, v0221, v1121); //62.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0122, v0221, v0212); //62.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0122, v1112, v1121); //62.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v0122, v1112, v0212); //62.3.3
        };
      };
    };        
  }
  else if (sB<(ValueType)1.0) //63
  {
    if(axis=='0') //63.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v2012, v2111, v2021); //63.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v2012, v2111, v1112); //63.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v2012, v1022, v2021); //63.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1121, v2012, v1022, v1112); //63.1.2
        };
      };
    }
    else if(axis=='1') //63.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2021, v1112, v2111, v1121); //63.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2021, v1112, v2111, v2012); //63.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2021, v1112, v1022, v1121); //63.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2021, v1112, v1022, v2012); //63.2.4
        };
      };
    }
    else //63.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1022, v1121, v2021); //63.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1022, v1121, v1112); //63.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1022, v2012, v2021); //63.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1022, v2012, v1112); //63.3.3
        };
      };
    };        
  }
  else if (sC<(ValueType)1.0) //64
  {
    if(axis=='0') //64.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v2102, v2201, v2111); //64.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v2102, v2201, v1202); //64.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v2102, v1112, v2111); //64.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1211, v2102, v1112, v1202); //64.1.2
        };
      };
    }
    else if(axis=='1') //64.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1202, v2201, v1211); //64.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1202, v2201, v2102); //64.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1202, v1112, v1211); //64.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1202, v1112, v2102); //64.2.4
        };
      };
    }
    else //64.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2201, v1112, v1211, v2111); //64.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2201, v1112, v1211, v1202); //64.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2201, v1112, v2102, v2111); //64.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2201, v1112, v2102, v1202); //64.3.3
        };
      };
    };        
  }
  else if (sD<(ValueType)1.0) //65
  {
    if(axis=='0') //65.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v2111, v2210, v2120); //65.1.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v2111, v2210, v1211); //65.1.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v2111, v1121, v2120); //65.1.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v1220, v2111, v1121, v1211); //65.1.2
        };
      };
    }
    else if(axis=='1') //65.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1211, v2210, v1220); //65.2.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1211, v2210, v2111); //65.2.3
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1211, v1121, v1220); //65.2.2
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2120, v1211, v1121, v2111); //65.2.4
        };
      };
    }
    else //65.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1121, v1220, v2120); //65.3.1
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1121, v1220, v1211); //65.3.2
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1121, v2111, v2120); //65.3.4
        }
        else
        {
          interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2210, v1121, v2111, v1211); //65.3.3
        };
      };
    };        
  }

  else //35
  {
    interpolationMinMax_3D<ValueType>(&Vmin, &Vmax, v2111, v1211, v1121, v1112);
  };

  if (result<Vmin)
  {
    result=Vmin;
  }
  else if (result>Vmax)
  {
    result=Vmax;
  };

  return(result);
}

//////////////////////ELSE////////////////////////////////////////////////////////////////////////


template <typename ValueType>
void interpolationMinMax_2D(ValueType * Min, ValueType * Max, 
const ValueType A, const ValueType B, const ValueType C)
{
  if (A<B)
  {
    (*Min)=A;
    (*Max)=B;
  }
  else
  {
    (*Min)=B;
    (*Max)=A;
  };
  if (C<(*Min))
  {
    (*Min)=C;
  }
  else if (C>(*Max))
  {
    (*Max)=C;
  };
  return;
}

template <typename ValueType>
void interpolationMinMax_3D(ValueType * Min, ValueType * Max, 
const ValueType A, const ValueType B, const ValueType C, const ValueType D)
{
  if (A<B)
  {
    (*Min)=A;
    (*Max)=B;
  }
  else
  {
    (*Min)=B;
    (*Max)=A;
  };
  if (C<(*Min))
  {
    (*Min)=C;
  }
  else if (C>(*Max))
  {
    (*Max)=C;
  };
  if (D<(*Min))
  {
    (*Min)=D;
  }
  else if (D>(*Max))
  {
    (*Max)=D;
  };
  return;
}

////////////////Координаты вершин маленького тетраэдра////////////////

template <typename ValueType>
void Nodes_Coord_3D(
ValueType * Coord,
const int N, 
const char Axis,
const ValueType * R) 
{
  int i;
  for(i=0; (i<12); i++)
  {
    Coord[i]=R[i+3];
  };
  if(N==2)
  {
    Nodes_Coord_3D_2<ValueType>(Coord, Axis,
                                             R[ 0], R[ 1], R[ 2],
                                             R[ 3], R[ 4], R[ 5],
                                             R[ 6], R[ 7], R[ 8],
                                             R[ 9], R[10], R[11],
                                             R[12], R[13], R[14]);
  }
  else if(N==3)
  {
    Nodes_Coord_3D_3<ValueType>(Coord, Axis, 
                                            R[ 0], R[ 1], R[ 2],
                                            R[ 3], R[ 4], R[ 5],
                                            R[ 6], R[ 7], R[ 8],
                                            R[ 9], R[10], R[11],
                                            R[12], R[13], R[14]);
  }
  else if(N==4)
  {
    Nodes_Coord_3D_4<ValueType>(Coord, Axis, 
                                             R[ 0], R[ 1], R[ 2],
                                             R[ 3], R[ 4], R[ 5],
                                             R[ 6], R[ 7], R[ 8],
                                             R[ 9], R[10], R[11],
                                             R[12], R[13], R[14]);
  }
  else if(N==5)
  {
    Nodes_Coord_3D_5<ValueType>(Coord, Axis, 
                                             R[ 0], R[ 1], R[ 2],
                                             R[ 3], R[ 4], R[ 5],
                                             R[ 6], R[ 7], R[ 8],
                                             R[ 9], R[10], R[11],
                                             R[12], R[13], R[14]);
  };
  return;
}

template <typename ValueType>
void Nodes_Coord_3D_2(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD)
{
  int i;
  int a [4];
  int b [4];
  int c [4];
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  sA=(ValueType)2.0*sA;
  sB=(ValueType)2.0*sB;
  sC=(ValueType)2.0*sC;
  sD=(ValueType)2.0*sD;  

  if (sA>(ValueType)1) //1
  {
    // v2000, v1100, v1010, v1001

    a[0]=2;
    b[0]=0;
    c[0]=0;

    a[1]=1;
    b[1]=1;
    c[1]=0;

    a[2]=1;
    b[2]=0;
    c[2]=1;

    a[3]=1;
    b[3]=0;
    c[3]=0;
  }
  else if (sB>(ValueType)1) //2
  {
    // v1100, v0200, v0110, v0101

    a[0]=1;
    b[0]=1;
    c[0]=0;

    a[1]=0;
    b[1]=2;
    c[1]=0;

    a[2]=0;
    b[2]=1;
    c[2]=1;

    a[3]=0;
    b[3]=1;
    c[3]=0;
  }
  else if (sC>(ValueType)1) //3
  {
    // v1010, v0110, v0020, v0011

    a[0]=1;
    b[0]=0;
    c[0]=1;

    a[1]=0;
    b[1]=1;
    c[1]=1;

    a[2]=0;
    b[2]=0;
    c[2]=2;

    a[3]=0;
    b[3]=0;
    c[3]=1;
  }
  else if (sD>(ValueType)1) //4
  {
    // v1001, v0101, v0011, v0002

    a[0]=1;
    b[0]=0;
    c[0]=0;

    a[1]=0;
    b[1]=1;
    c[1]=0;

    a[2]=0;
    b[2]=0;
    c[2]=1;

    a[3]=0;
    b[3]=0;
    c[3]=0;
  }
  else //5
  {
    if(axis=='0') //5.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0110, v1001, v1100, v1010 - 5.1.1

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0110, v1001, v1100, v0101 - 5.1.3

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0110, v1001, v0011, v1010 - 5.1.4

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0110, v1001, v0011, v0101 - 5.1.2

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //5.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1010, v0101, v1100, v0110 - 5.2.1

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1010, v0101, v1100, v1001 - 5.2.3

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1010, v0101, v0011, v0110 - 5.2.2

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1010, v0101, v0011, v1001 - 5.2.4

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //5.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1100, v0011, v0110, v1010 - 5.3.1

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1100, v0011, v0110, v0101 - //5.3.2

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1100, v0011, v1001, v1010) - 5.3.4

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1100, v0011, v1001, v0101 - 5.3.3

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  };

  for(i=0; (i<4); i++)
  {
    Coord[i*3  ]=xD+(xA-xD)*a[i]+(xB-xD)*b[i]+(xC-xD)*c[i];
    Coord[i*3+1]=yD+(yA-yD)*a[i]+(yB-yD)*b[i]+(yC-yD)*c[i];
    Coord[i*3+2]=zD+(zA-zD)*a[i]+(zB-zD)*b[i]+(zC-zD)*c[i];
  };

  return;
}

template <typename ValueType>
void Nodes_Coord_3D_3(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD)
{
  int i;
  int a [4];
  int b [4];
  int c [4];
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  sA=(ValueType)3.0*sA;
  sB=(ValueType)3.0*sB;
  sC=(ValueType)3.0*sC;
  sD=(ValueType)3.0*sD;  

  if (sA>(ValueType)2) //1
  {
    // v3000, v2100, v2010, v2001

    a[0]=3;
    b[0]=0;
    c[0]=0;

    a[1]=2;
    b[1]=1;
    c[1]=0;

    a[2]=2;
    b[2]=0;
    c[2]=1;

    a[3]=2;
    b[3]=0;
    c[3]=0;
  }
  else if(sB>(ValueType)2) //2
  {
    // v1200, v0300, v0210, v0201

    a[0]=1;
    b[0]=2;
    c[0]=0;

    a[1]=0;
    b[1]=3;
    c[1]=0;

    a[2]=0;
    b[2]=2;
    c[2]=1;

    a[3]=0;
    b[3]=2;
    c[3]=0;
  }
  else if(sC>(ValueType)2) //3
  {
    // v1020, v0120, v0030, v0021

    a[0]=1;
    b[0]=0;
    c[0]=2;

    a[1]=0;
    b[1]=1;
    c[1]=2;

    a[2]=0;
    b[2]=0;
    c[2]=3;

    a[3]=0;
    b[3]=0;
    c[3]=2;
  }
  else if(sD>(ValueType)2) //4
  {
    // v1002, v0102, v0012, v0003

    a[0]=1;
    b[0]=0;
    c[0]=0;

    a[1]=0;
    b[1]=1;
    c[1]=0;

    a[2]=0;
    b[2]=0;
    c[2]=1;

    a[3]=0;
    b[3]=0;
    c[3]=0;
  }
  else if((sA>(ValueType)1)&&(sB>(ValueType)1)) //5
  {
    // v2100, v1200, v1110, v1101

    a[0]=2;
    b[0]=1;
    c[0]=0;

    a[1]=1;
    b[1]=2;
    c[1]=0;

    a[2]=1;
    b[2]=1;
    c[2]=1;

    a[3]=1;
    b[3]=1;
    c[3]=0;
  }
  else if((sB>(ValueType)1)&&(sC>(ValueType)1)) //6
  {
    // v1110, v0210, v0120, v0111 

    a[0]=1;
    b[0]=1;
    c[0]=1;

    a[1]=0;
    b[1]=2;
    c[1]=1;

    a[2]=0;
    b[2]=1;
    c[2]=2;

    a[3]=0;
    b[3]=1;
    c[3]=1;
  }
  else if((sC>(ValueType)1)&&(sD>(ValueType)1)) //7
  {
    // v1011, v0111, v0021, v0012

    a[0]=1;
    b[0]=0;
    c[0]=1;

    a[1]=0;
    b[1]=1;
    c[1]=1;

    a[2]=0;
    b[2]=0;
    c[2]=2;

    a[3]=0;
    b[3]=0;
    c[3]=1;
  }
  else if((sA>(ValueType)1)&&(sD>(ValueType)1)) //8
  {
    // v2001, v1101, v1011, v1002 

    a[0]=2;
    b[0]=0;
    c[0]=0;

    a[1]=1;
    b[1]=1;
    c[1]=0;

    a[2]=1;
    b[2]=0;
    c[2]=1;

    a[3]=1;
    b[3]=0;
    c[3]=0;
  }
  else if((sA>(ValueType)1)&&(sC>(ValueType)1)) //9
  {
    // v2010, v1110, v1020, v1011 

    a[0]=2;
    b[0]=0;
    c[0]=1;

    a[1]=1;
    b[1]=1;
    c[1]=1;

    a[2]=1;
    b[2]=0;
    c[2]=2;

    a[3]=1;
    b[3]=0;
    c[3]=1;
  }
  else if((sB>(ValueType)1)&&(sD>(ValueType)1)) //10
  {
    // v1101, v0201, v0111, v0102

    a[0]=1;
    b[0]=1;
    c[0]=0;

    a[1]=0;
    b[1]=2;
    c[1]=0;

    a[2]=0;
    b[2]=1;
    c[2]=1;

    a[3]=0;
    b[3]=1;
    c[3]=0;
  }
  else if(sA>(ValueType)1) //11
  {
    if(axis=='0') //11.1
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1110, v2001, v2100, v2010 - 11.1.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1110, v2001, v2100, v1101 - 11.1.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1110, v2001, v1011, v2010 - 11.1.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1110, v2001, v1011, v1101 - 11.1.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //11.2
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2010, v1101, v2100, v1110 - 11.2.1

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2010, v1101, v2100, v2001 - 11.2.3

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2010, v1101, v1011, v1110 - 11.2.2

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2010, v1101, v1011, v2001 - 11.2.4

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //11.3
    {
      s1=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2100, v1011, v1110, v2010) - 11.3.1

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2100, v1011, v1110, v1101 - 11.3.2

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2100, v1011, v2001, v2010 - 11.3.4

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2100, v1011, v2001, v1101 - 11.3.3

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    };     
  }
  else if(sB>(ValueType)1) //12
  {
    if(axis=='0') //12.1
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0210, v1101, v1200, v1110 - 12.1.1

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v0210, v1101, v1200, v0201 - 12.1.3

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0210, v1101, v0111, v1110 - 12.1.4

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v0210, v1101, v0111, v0201 - 12.1.2

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //12.2
    {
      s1=(sA+sB-sC-sD-(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1110, v0201, v1200, v0210 - 12.2.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1110, v0201, v1200, v1101 - 12.2.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1110, v0201, v0111, v0210 - 12.2.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1110, v0201, v0111, v1101 - 12.2.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else //12.3
    {
      s1=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1200, v0111, v0210, v1110 - 12.3.1

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1200, v0111, v0210, v0201 - 12.3.2

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1200, v0111, v1101, v1110 - 12.3.4

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1200, v0111, v1101, v0201 - 12.3.3

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      };
    };  
  }
  else if(sC>(ValueType)1) //13
  {
    if(axis=='0') //13.1
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0120, v1011, v1110, v1020 - 13.1.1

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v0120, v1011, v1110, v0111 - 13.1.3

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0120, v1011, v0021, v1020 - 13.1.4

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v0120, v1011, v0021, v0111 - 13.1.2

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //13.2
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1020, v0111, v1110, v0120 - 13.2.1

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1020, v0111, v1110, v1011 - 13.2.3

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1020, v0111, v0021, v0120 - 13.2.2

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1020, v0111, v0021, v1011); //13.2.4

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        };
      };
    }
    else //13.3
    {
      s1=(sB+sC-sD-sA-(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD-(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1110, v0021, v0120, v1020 - 13.3.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1110, v0021, v0120, v0111 - 13.3.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1110, v0021, v1011, v1020 - 13.3.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1110, v0021, v1011, v0111 - 13.3.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      };
    };  
  }      
  else if(sD>(ValueType)1)//14
  {
    if(axis=='0') //14.1
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0111, v1002, v1101, v1011 - 14.1.1

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0111, v1002, v1101, v0102 - 14.1.3

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0111, v1002, v0012, v1011 - 14.1.4

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0111, v1002, v0012, v0102 - 14.1.2

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //14.2
    {
      s1=(sA+sB-sC-sD+(ValueType)1.0)/(ValueType)2.0;
      s2=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1011, v0102, v1101, v0111 - 14.2.1

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1011, v0102, v1101, v1002 - 14.2.3

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1011, v0102, v0012, v0111 - 14.2.2

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1011, v0102, v0012, v1002 - 14.2.4

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //14.3
    {
      s1=(sB+sC-sD-sA+(ValueType)1.0)/(ValueType)2.0;
      s2=(sA+sC-sB-sD+(ValueType)1.0)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1101, v0012, v0111, v1011 - 14.3.1

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1101, v0012, v0111, v0102 - 14.3.2

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1101, v0012, v1002, v1011 - 14.3.4

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1101, v0012, v1002, v0102 - 14.3.3

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      };
    };
  }  
  else //15
  {
    // v0111, v1011, v1101, v1110

    a[0]=0;
    b[0]=1;
    c[0]=1;

    a[1]=1;
    b[1]=0;
    c[1]=1;

    a[2]=1;
    b[2]=1;
    c[2]=0;

    a[3]=1;
    b[3]=1;
    c[3]=1;
  };

  for(i=0; (i<4); i++)
  {
    Coord[i*3  ]=xD+(xA-xD)*a[i]+(xB-xD)*b[i]+(xC-xD)*c[i];
    Coord[i*3+1]=yD+(yA-yD)*a[i]+(yB-yD)*b[i]+(yC-yD)*c[i];
    Coord[i*3+2]=zD+(zA-zD)*a[i]+(zB-zD)*b[i]+(zC-zD)*c[i];
  };

  return;
}

template <typename ValueType>
void Nodes_Coord_3D_4(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD)
{
  int i;
  int a [4];
  int b [4];
  int c [4];
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  sA=(ValueType)4.0*sA;
  sB=(ValueType)4.0*sB;
  sC=(ValueType)4.0*sC;
  sD=(ValueType)4.0*sD;  

  if (sA>(ValueType)3.0) //1
  {
    // v4000, v3100, v3010, v3001

    a[0]=4;
    b[0]=0;
    c[0]=0;

    a[1]=3;
    b[1]=1;
    c[1]=0;

    a[2]=3;
    b[2]=0;
    c[2]=1;

    a[3]=3;
    b[3]=0;
    c[3]=0;
  }
  else if (sB>(ValueType)3.0) //2
  {
    // v1300, v0400, v0310, v0301

    a[0]=1;
    b[0]=3;
    c[0]=0;

    a[1]=0;
    b[1]=4;
    c[1]=0;

    a[2]=0;
    b[2]=3;
    c[2]=1;

    a[3]=0;
    b[3]=3;
    c[3]=0;
  }
  else if (sC>(ValueType)3.0) //3
  {
    // v1030, v0130, v0040, v0031

    a[0]=1;
    b[0]=0;
    c[0]=3;

    a[1]=0;
    b[1]=1;
    c[1]=3;

    a[2]=0;
    b[2]=0;
    c[2]=4;

    a[3]=0;
    b[3]=0;
    c[3]=3;
  }
  else if (sD>(ValueType)3.0) //4
  {
    // v1003, v0103, v0013, v0004

    a[0]=1;
    b[0]=0;
    c[0]=0;

    a[1]=0;
    b[1]=1;
    c[1]=0;

    a[2]=0;
    b[2]=0;
    c[2]=1;

    a[3]=0;
    b[3]=0;
    c[3]=0;
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)) //9
  {
    // v3100, v2200, v2110, v2101

    a[0]=3;
    b[0]=1;
    c[0]=0;

    a[1]=2;
    b[1]=2;
    c[1]=0;

    a[2]=2;
    b[2]=1;
    c[2]=1;

    a[3]=2;
    b[3]=1;
    c[3]=0;
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //10
  {
    // v1210, v0310, v0220, v0211

    a[0]=1;
    b[0]=2;
    c[0]=1;

    a[1]=0;
    b[1]=3;
    c[1]=1;

    a[2]=0;
    b[2]=2;
    c[2]=2;

    a[3]=0;
    b[3]=2;
    c[3]=1;
  }
  else if ((sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //11
  {
    // v1021, v0121, v0031, v0022

    a[0]=1;
    b[0]=0;
    c[0]=2;

    a[1]=0;
    b[1]=1;
    c[1]=2;

    a[2]=0;
    b[2]=0;
    c[2]=3;

    a[3]=0;
    b[3]=0;
    c[3]=2;
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //12
  {
    // v2002, v1102, v1012, v1003

    a[0]=2;
    b[0]=0;
    c[0]=0;

    a[1]=1;
    b[1]=1;
    c[1]=0;

    a[2]=1;
    b[2]=0;
    c[2]=1;

    a[3]=1;
    b[3]=0;
    c[3]=0;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)) //13
  {
    // v2200, v1300, v1210, v1201

    a[0]=2;
    b[0]=2;
    c[0]=0;

    a[1]=1;
    b[1]=3;
    c[1]=0;

    a[2]=1;
    b[2]=2;
    c[2]=1;

    a[3]=1;
    b[3]=2;
    c[3]=0;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //14
  {
    // v1120, v0220, v0130, v0121

    a[0]=1;
    b[0]=1;
    c[0]=2;

    a[1]=0;
    b[1]=2;
    c[1]=2;

    a[2]=0;
    b[2]=1;
    c[2]=3;

    a[3]=0;
    b[3]=1;
    c[3]=2;
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //15
  {
    // v1012, v0112, v0022, v0013

    a[0]=1;
    b[0]=0;
    c[0]=1;

    a[1]=0;
    b[1]=1;
    c[1]=1;

    a[2]=0;
    b[2]=0;
    c[2]=2;

    a[3]=0;
    b[3]=0;
    c[3]=1;
  }
  else if ((sA>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //16
  {
    // v3001, v2101, v2011, v2002

    a[0]=3;
    b[0]=0;
    c[0]=0;

    a[1]=2;
    b[1]=1;
    c[1]=0;

    a[2]=2;
    b[2]=0;
    c[2]=1;

    a[3]=2;
    b[3]=0;
    c[3]=0;
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //17
  {
    // v3010, v2110, v2020, v2011

    a[0]=3;
    b[0]=0;
    c[0]=1;

    a[1]=2;
    b[1]=1;
    c[1]=1;

    a[2]=2;
    b[2]=0;
    c[2]=2;

    a[3]=2;
    b[3]=0;
    c[3]=1;
  }
  else if ((sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //18
  {
    // v1201, v0301, v0211, v0202

    a[0]=1;
    b[0]=2;
    c[0]=0;

    a[1]=0;
    b[1]=3;
    c[1]=0;

    a[2]=0;
    b[2]=2;
    c[2]=1;

    a[3]=0;
    b[3]=2;
    c[3]=0;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //19
  {
    // v2020, v1120, v1030, v1021

    a[0]=2;
    b[0]=0;
    c[0]=2;

    a[1]=1;
    b[1]=1;
    c[1]=2;

    a[2]=1;
    b[2]=0;
    c[2]=3;

    a[3]=1;
    b[3]=0;
    c[3]=2;
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //20
  {
    // v1102, v0202, v0112, v0103

    a[0]=1;
    b[0]=1;
    c[0]=0;

    a[1]=0;
    b[1]=2;
    c[1]=0;

    a[2]=0;
    b[2]=1;
    c[2]=1;

    a[3]=0;
    b[3]=1;
    c[3]=0;
  }

  else if (sA>(ValueType)2.0) //25
  {
    if(axis=='0') //25.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2110, v3001, v3100, v3010 - 25.1.1

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2110, v3001, v3100, v2101 - 25.1.3

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2110, v3001, v2011, v3010 - 25.1.4

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2110, v3001, v2011, v2101 - 25.1.2

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //25.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3010, v2101, v3100, v2110 - 25.2.1

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v3010, v2101, v3100, v3001 - 25.2.3

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=3;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3010, v2101, v2011, v2110 - 25.2.2

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v3010, v2101, v2011, v3001 - 25.2.4

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //25.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3100, v2011, v2110, v3010 - 25.3.1

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v3100, v2011, v2110, v2101 - 25.3.2

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3100, v2011, v3001, v3010 - 25.3.4

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=0;
          c[2]=0;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v3100, v2011, v3001, v2101 - 25.3.3

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=0;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  }
  else if (sB>(ValueType)2.0) //26
  {
    if(axis=='0') //26.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0310, v1201, v1300, v1210 - 26.1.1

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v0310, v1201, v1300, v0301 - 26.1.3

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0310, v1201, v0211, v1210 - 26.1.4

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v0310, v1201, v0211, v0301 - 26.1.2

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //26.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1210, v0301, v1300, v0310 - 26.2.1

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=0;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v1210, v0301, v1300, v1201 - 26.2.3 

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1210, v0301, v0211, v0310 - 26.2.2

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v1210, v0301, v0211, v1201 - 26.2.4

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else //26.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1300, v0211, v0310, v1210 - 26.3.1

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1300, v0211, v0310, v0301 - 26.3.2

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1300, v0211, v1201, v1210 - 26.3.4

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1300, v0211, v1201, v0301 - 26.3.3

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      };
    };        
  }
  else if (sC>(ValueType)2.0) //27
  {
    if(axis=='0') //27.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0130, v1021, v1120, v1030 - 27.1.1

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v0130, v1021, v1120, v0121 - 27.1.3

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0130, v1021, v0031, v1030 - 27.1.4

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v0130, v1021, v0031, v0121 - 27.1.2

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      };
    }
    else if(axis=='1') //27.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1030, v0121, v1120, v0130 - 27.2.1
 
          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=3;
       }
        else
        {
          // v1030, v0121, v1120, v1021 - 27.2.3

          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1030, v0121, v0031, v0130 - 27.2.2

          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v1030, v0121, v0031, v1021 - 27.2.4

          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        };
      };
    }
    else //27.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1120, v0031, v0130, v1030 - 27.3.1

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v1120, v0031, v0130, v0121 - 27.3.2

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1120, v0031, v1021, v1030 - 27.3.4

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v1120, v0031, v1021, v0121 - 27.3.3

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      };
    };        
  }
  else if (sD>(ValueType)2.0) //28
  {
    if(axis=='0') //28.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0112, v1003, v1102, v1012 - 28.1.1

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0112, v1003, v1102, v0103 - 28.1.3

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0112, v1003, v0013, v1012 - 28.1.4

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0112, v1003, v0013, v0103 - 28.1.2
 
          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
       };
      };
    }
    else if(axis=='1') //28.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1012, v0103, v1102, v0112 - 28.2.1

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1012, v0103, v1102, v1003 - 28.2.3

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1012, v0103, v0013, v0112 - 28.2.2

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1012, v0103, v0013, v1003 - 28.2.4

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //28.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1102, v0013, v0112, v1012 - 28.3.1

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1102, v0013, v0112, v0103 - 28.3.2

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1102, v0013, v1003, v1012 - 28.3.4

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1102, v0013, v1003, v0103 - 28.3.3

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  }

  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //21
  {
    // v1111, v2011, v2101, v2110 

    a[0]=1;
    b[0]=1;
    c[0]=1;

    a[1]=2;
    b[1]=0;
    c[1]=1;

    a[2]=2;
    b[2]=1;
    c[2]=0;

    a[3]=2;
    b[3]=1;
    c[3]=1;
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //22
  {
    // v0211, v1111, v1201, v1210 

    a[0]=0;
    b[0]=2;
    c[0]=1;

    a[1]=1;
    b[1]=1;
    c[1]=1;

    a[2]=1;
    b[2]=2;
    c[2]=0;

    a[3]=1;
    b[3]=2;
    c[3]=1;
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //23
  {
    // v0121, v1021, v1111, v1120 

    a[0]=0;
    b[0]=1;
    c[0]=2;

    a[1]=1;
    b[1]=0;
    c[1]=2;

    a[2]=1;
    b[2]=1;
    c[2]=1;

    a[3]=1;
    b[3]=1;
    c[3]=2;
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //24
  {
    // v0112, v1012, v1102, v1111 

    a[0]=0;
    b[0]=1;
    c[0]=1;

    a[1]=1;
    b[1]=0;
    c[1]=1;

    a[2]=1;
    b[2]=1;
    c[2]=0;

    a[3]=1;
    b[3]=1;
    c[3]=1;
  }

  else if ((sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //29
  {
    if(axis=='0') //29.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1210, v2101, v2200, v2110 - 29.1.1

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1210, v2101, v2200, v1201 - 29.1.3

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1210, v2101, v1111, v2110 - 29.1.4

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1210, v2101, v1111, v1201 - 29.1.2

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //29.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2110, v1201, v2200, v1210 - 29.2.1

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v2110, v1201, v2200, v2101 - 29.2.3

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2110, v1201, v1111, v1210 - 29.2.2

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v2110, v1201, v1111, v2101 - 29.2.4

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else //29.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2200, v1111, v1210, v2110 - 29.3.1

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2200, v1111, v1210, v1201 - 29.3.2

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2200, v1111, v2101, v2110 - 29.3.4

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2200, v1111, v2101, v1201); //29.3.3

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sD<(ValueType)1.0)) //30
  {
    if(axis=='0') //30.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0220, v1111, v1210, v1120 - 30.1.1

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v0220, v1111, v1210, v0211 - 30.1.3

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0220, v1111, v0121, v1120 - 30.1.4

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v0220, v1111, v0121, v0211 - 30.1.2

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //30.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1120, v0211, v1210, v0220 - 30.2.1

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v1120, v0211, v1210, v1111 - 30.2.3

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1120, v0211, v0121, v0220 - 30.2.2

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v1120, v0211, v0121, v1111 - 30.2.4

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else //30.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1210, v0121, v0220, v1120 - 30.3.1

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1210, v0121, v0220, v0211 - 30.3.2

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1210, v0121, v1111, v1120 - 30.3.4

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1210, v0121, v1111, v0211 - 30.3.3

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //31
  {
    if(axis=='0') //31.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0121, v1012, v1111, v1021 - 31.1.1

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v0121, v1012, v1111, v0112 - 31.1.3

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0121, v1012, v0022, v1021 - 31.1.4

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v0121, v1012, v0022, v0112 - 31.1.2

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //31.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1021, v0112, v1111, v0121 - 31.2.1

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1021, v0112, v1111, v1012 - 31.2.3

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1021, v0112, v0022, v0121 - 31.2.2

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1021, v0112, v0022, v1012 - 31.2.4

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        };
      };
    }
    else //31.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1111, v0022, v0121, v1021 - 31.3.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1111, v0022, v0121, v0112 - 31.3.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1111, v0022, v1012, v1021 - 31.3.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1111, v0022, v1012, v0112 - 31.3.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      };
    };        
  }
  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //32
  {
    if(axis=='0') //32.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1111, v2002, v2101, v2011 - 32.1.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1111, v2002, v2101, v1102 - 32.1.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1111, v2002, v1012, v2011 - 32.1.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1111, v2002, v1012, v1102); //32.1.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //32.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2011, v1102, v2101, v1111 - 32.2.1

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2011, v1102, v2101, v2002); //32.2.3

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2011, v1102, v1012, v1111 - 32.2.2

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2011, v1102, v1012, v2002); //32.2.4

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //32.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2101, v1012, v1111, v2011 - 32.3.1

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2101, v1012, v1111, v1102 - 32.3.2

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2101, v1012, v2002, v2011 - 32.3.4

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2101, v1012, v2002, v1102 - 32.3.3

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  }
  else if ((sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //33
  {
    if(axis=='0') //33.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1120, v2011, v2200, v2110 - 33.1.1

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1120, v2011, v2200, v1201 - 33.1.3

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1120, v2011, v1111, v2110 - 33.1.4

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1120, v2011, v1111, v1201 - 33.1.2

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //33.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2020, v1111, v2110, v1120 - 33.2.1

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v2020, v1111, v2110, v2011 - 33.2.3

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2020, v1111, v1021, v1120 - 33.2.2

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v2020, v1111, v1021, v2011 - 33.2.4

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        };
      };
    }
    else //33.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2110, v1021, v1120, v2020 - 33.3.1

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v2110, v1021, v1120, v1111 - 33.3.2

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2110, v1021, v2011, v2020 - 33.3.4

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v2110, v1021, v2011, v1111 - 33.3.3

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      };
    };        
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //34
  {
    if(axis=='0') //34.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0211, v1102, v1201, v1111 - 34.1.1

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v0211, v1102, v1201, v0202 - 34.1.3

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0211, v1102, v0112, v1111 - 34.1.4

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v0211, v1102, v0112, v0202 - 34.1.2

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //34.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0;
      s2=(sB+sC-sD-sA)/(ValueType)2.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1111, v0202, v1201, v0211 - 34.2.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1111, v0202, v1201, v1102 - 34.2.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1111, v0202, v0112, v0211); //34.2.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1111, v0202, v0112, v1102 - 34.2.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else //34.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.0;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1201, v0112, v0211, v1111 - 34.3.1

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1201, v0112, v0211, v0202 - 34.3.2

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1201, v0112, v1102, v1111 - 34.3.4

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1201, v0112, v1102, v0202 - 34.3.3

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      };
    };        
  }

  else if (sA<(ValueType)1.0) //5
  {
    // v1111, v0211, v0121, v0112 

    a[0]=1;
    b[0]=1;
    c[0]=1;

    a[1]=0;
    b[1]=2;
    c[1]=1;

    a[2]=0;
    b[2]=1;
    c[2]=2;

    a[3]=0;
    b[3]=1;
    c[3]=1;
  }
  else if (sB<(ValueType)1.0) //6
  {
    // v2011, v1111, v1021, v1012 

    a[0]=2;
    b[0]=0;
    c[0]=1;

    a[1]=1;
    b[1]=1;
    c[1]=1;

    a[2]=1;
    b[2]=0;
    c[2]=2;

    a[3]=1;
    b[3]=0;
    c[3]=1;
  }
  else if (sC<(ValueType)1.0) //7
  {
    // v2101, v1201, v1111, v1102 

    a[0]=2;
    b[0]=1;
    c[0]=0;

    a[1]=1;
    b[1]=2;
    c[1]=0;

    a[2]=1;
    b[2]=1;
    c[2]=1;

    a[3]=1;
    b[3]=1;
    c[3]=0;
  }
  else //8
  {
    // v2110, v1210, v1120, v1111 

    a[0]=2;
    b[0]=1;
    c[0]=1;

    a[1]=1;
    b[1]=2;
    c[1]=1;

    a[2]=1;
    b[2]=1;
    c[2]=2;

    a[3]=1;
    b[3]=1;
    c[3]=1;
  };

  for(i=0; (i<4); i++)
  {
    Coord[i*3  ]=xD+(xA-xD)*a[i]+(xB-xD)*b[i]+(xC-xD)*c[i];
    Coord[i*3+1]=yD+(yA-yD)*a[i]+(yB-yD)*b[i]+(yC-yD)*c[i];
    Coord[i*3+2]=zD+(zA-zD)*a[i]+(zB-zD)*b[i]+(zC-zD)*c[i];
  };

  return;
}


template <typename ValueType>
void Nodes_Coord_3D_5(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD)
{
  int i;
  int a [4];
  int b [4];
  int c [4];
  ValueType sA, sB, sC, sD;
  ValueType s1, s2;

  interpolationInitial_3D<ValueType>(x, y, z, 
                                     xA, yA, zA, 
                                     xB, yB, zB, 
                                     xC, yC, zC, 
                                     xD, yD, zD, 
                                     &sA, &sB, &sC, &sD);

  sA=(ValueType)5.0*sA;
  sB=(ValueType)5.0*sB;
  sC=(ValueType)5.0*sC;
  sD=(ValueType)5.0*sD;  

  if (sA>(ValueType)4.0) //1
  {
    // v5000, v4100, v4010, v4001 

    a[0]=5;
    b[0]=0;
    c[0]=0;

    a[1]=4;
    b[1]=1;
    c[1]=0;

    a[2]=4;
    b[2]=0;
    c[2]=1;

    a[3]=4;
    b[3]=0;
    c[3]=0;
  }
  else if (sB>(ValueType)4.0) //2
  {
    // v1400, v0500, v0410, v0401 

    a[0]=1;
    b[0]=4;
    c[0]=0;

    a[1]=0;
    b[1]=5;
    c[1]=0;

    a[2]=0;
    b[2]=4;
    c[2]=1;

    a[3]=0;
    b[3]=4;
    c[3]=0;
  }
  else if (sC>(ValueType)4.0) //3
  {
    // v1040, v0140, v0050, v0041 

    a[0]=1;
    b[0]=0;
    c[0]=4;

    a[1]=0;
    b[1]=1;
    c[1]=4;

    a[2]=0;
    b[2]=0;
    c[2]=5;

    a[3]=0;
    b[3]=0;
    c[3]=4;
  }
  else if (sD>(ValueType)4.0) //4
  {
    // v1004, v0104, v0014, v0005 

    a[0]=1;
    b[0]=0;
    c[0]=0;

    a[1]=0;
    b[1]=1;
    c[1]=0;

    a[2]=0;
    b[2]=0;
    c[2]=1;

    a[3]=0;
    b[3]=0;
    c[3]=0;
  }

  else if ((sA>(ValueType)3.0)&&(sB>=(ValueType)1.0)) //11
  {
    // v4100, v3200, v3110, v3101 

    a[0]=4;
    b[0]=1;
    c[0]=0;

    a[1]=3;
    b[1]=2;
    c[1]=0;

    a[2]=3;
    b[2]=1;
    c[2]=1;

    a[3]=3;
    b[3]=1;
    c[3]=0;
  }
  else if ((sB>(ValueType)3.0)&&(sC>=(ValueType)1.0)) //12
  {
    // v1310, v0410, v0320, v0311 

    a[0]=1;
    b[0]=3;
    c[0]=1;

    a[1]=0;
    b[1]=4;
    c[1]=1;

    a[2]=0;
    b[2]=3;
    c[2]=2;

    a[3]=0;
    b[3]=3;
    c[3]=1;
  }
  else if ((sC>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //13
  {
    // v1031, v0131, v0041, v0032 

    a[0]=1;
    b[0]=0;
    c[0]=3;

    a[1]=0;
    b[1]=1;
    c[1]=3;

    a[2]=0;
    b[2]=0;
    c[2]=4;

    a[3]=0;
    b[3]=0;
    c[3]=3;
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //14
  {
    // v2003, v1103, v1013, v1004 

    a[0]=2;
    b[0]=0;
    c[0]=0;

    a[1]=1;
    b[1]=1;
    c[1]=0;

    a[2]=1;
    b[2]=0;
    c[2]=1;

    a[3]=1;
    b[3]=0;
    c[3]=0;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)3.0)) //15
  {
    // v2300, v1400, v1310, v1301 

    a[0]=2;
    b[0]=3;
    c[0]=0;

    a[1]=1;
    b[1]=4;
    c[1]=0;

    a[2]=1;
    b[2]=3;
    c[2]=1;

    a[3]=1;
    b[3]=3;
    c[3]=0;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)3.0)) //16
  {
    // v1130, v0230, v0140, v0131 

    a[0]=1;
    b[0]=1;
    c[0]=3;

    a[1]=0;
    b[1]=2;
    c[1]=3;

    a[2]=0;
    b[2]=1;
    c[2]=4;

    a[3]=0;
    b[3]=1;
    c[3]=3;
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //17
  {
    // v1013, v0113, v0023, v0014 

    a[0]=1;
    b[0]=0;
    c[0]=1;

    a[1]=0;
    b[1]=1;
    c[1]=1;

    a[2]=0;
    b[2]=0;
    c[2]=2;

    a[3]=0;
    b[3]=0;
    c[3]=1;
  }
  else if ((sA>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //18
  {
    // v4001, v3101, v3011, v3002 

    a[0]=4;
    b[0]=0;
    c[0]=0;

    a[1]=3;
    b[1]=1;
    c[1]=0;

    a[2]=3;
    b[2]=0;
    c[2]=1;

    a[3]=3;
    b[3]=0;
    c[3]=0;
  }
  else if ((sA>(ValueType)3.0)&&(sC>=(ValueType)1.0)) //19
  {
    // v4010, v3110, v3020, v3011 

    a[0]=4;
    b[0]=0;
    c[0]=1;

    a[1]=3;
    b[1]=1;
    c[1]=1;

    a[2]=3;
    b[2]=0;
    c[2]=2;

    a[3]=3;
    b[3]=0;
    c[3]=1;
  }
  else if ((sB>(ValueType)3.0)&&(sD>=(ValueType)1.0)) //20
  {
    // v1301, v0401, v0311, v0302 

    a[0]=1;
    b[0]=3;
    c[0]=0;

    a[1]=0;
    b[1]=4;
    c[1]=0;

    a[2]=0;
    b[2]=3;
    c[2]=1;

    a[3]=0;
    b[3]=3;
    c[3]=0;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)3.0)) //21
  {
    // v2030, v1130, v1040, v1031 

    a[0]=2;
    b[0]=0;
    c[0]=3;

    a[1]=1;
    b[1]=1;
    c[1]=3;

    a[2]=1;
    b[2]=0;
    c[2]=4;

    a[3]=1;
    b[3]=0;
    c[3]=3;
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)3.0)) //22
  {
    // v1103, v0203, v0113, v0104 

    a[0]=1;
    b[0]=1;
    c[0]=0;

    a[1]=0;
    b[1]=2;
    c[1]=0;

    a[2]=0;
    b[2]=1;
    c[2]=1;

    a[3]=0;
    b[3]=1;
    c[3]=0;
  }

  else if (sA>(ValueType)3.0) //46
  {
    if(axis=='0') //46.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3110, v4001, v4100, v4010 - 46.1.1

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=4;
          b[1]=0;
          c[1]=0;

          a[2]=4;
          b[2]=1;
          c[2]=0;

          a[3]=4;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v3110, v4001, v4100, v3101 - 46.1.3

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=4;
          b[1]=0;
          c[1]=0;

          a[2]=4;
          b[2]=1;
          c[2]=0;

          a[3]=3;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3110, v4001, v3011, v4010 - 46.1.4

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=4;
          b[1]=0;
          c[1]=0;

          a[2]=3;
          b[2]=0;
          c[2]=1;

          a[3]=4;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v3110, v4001, v3011, v3101); //46.1.2

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=4;
          b[1]=0;
          c[1]=0;

          a[2]=3;
          b[2]=0;
          c[2]=1;

          a[3]=3;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //46.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v4010, v3101, v4100, v3110 - 46.2.1

          a[0]=4;
          b[0]=0;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=4;
          b[2]=1;
          c[2]=0;

          a[3]=3;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v4010, v3101, v4100, v4001 - 46.2.3

          a[0]=4;
          b[0]=0;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=4;
          b[2]=1;
          c[2]=0;

          a[3]=4;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v4010, v3101, v3011, v3110 - 46.2.2

          a[0]=4;
          b[0]=0;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=0;
          c[2]=1;

          a[3]=3;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v4010, v3101, v3011, v4001 - 46.2.4

          a[0]=4;
          b[0]=0;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=0;
          c[2]=1;

          a[3]=4;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //46.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v4100, v3011, v3110, v4010 - 46.3.1

          a[0]=4;
          b[0]=1;
          c[0]=0;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=1;

          a[3]=4;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v4100, v3011, v3110, v3101 - 46.3.2

          a[0]=4;
          b[0]=1;
          c[0]=0;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=1;

          a[3]=3;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v4100, v3011, v4001, v4010 - 46.3.4

          a[0]=4;
          b[0]=1;
          c[0]=0;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=4;
          b[2]=0;
          c[2]=0;

          a[3]=4;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v4100, v3011, v4001, v3101 - 46.3.3

          a[0]=4;
          b[0]=1;
          c[0]=0;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=4;
          b[2]=0;
          c[2]=0;

          a[3]=3;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  }
  else if (sB>(ValueType)3.0) //47
  {
    if(axis=='0') //47.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0410, v1301, v1400, v1310 - 47.1.1

          a[0]=0;
          b[0]=4;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=4;
          c[2]=0;

          a[3]=1;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v0410, v1301, v1400, v0401 - 47.1.3

          a[0]=0;
          b[0]=4;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=4;
          c[2]=0;

          a[3]=0;
          b[3]=4;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0410, v1301, v0311, v1310 - 47.1.4

          a[0]=0;
          b[0]=4;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=1;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v0410, v1301, v0311, v0401 - 47.1.2

          a[0]=0;
          b[0]=4;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=0;
          b[3]=4;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //47.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1310, v0401, v1400, v0410 - 47.2.1

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=4;
          c[1]=0;

          a[2]=1;
          b[2]=4;
          c[2]=0;

          a[3]=0;
          b[3]=4;
          c[3]=1;
        }
        else
        {
          // v1310, v0401, v1400, v1301 - 47.2.3

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=4;
          c[1]=0;

          a[2]=1;
          b[2]=4;
          c[2]=0;

          a[3]=1;
          b[3]=3;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1310, v0401, v0311, v0410 - 47.2.2

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=4;
          c[1]=0;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=0;
          b[3]=4;
          c[3]=1;
        }
        else
        {
          // v1310, v0401, v0311, v1301 - 47.2.4

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=4;
          c[1]=0;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=1;
          b[3]=3;
          c[3]=0;
        };
      };
    }
    else //47.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1400, v0311, v0410, v1310 - 47.3.1

          a[0]=1;
          b[0]=4;
          c[0]=0;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=0;
          b[2]=4;
          c[2]=1;

          a[3]=1;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v1400, v0311, v0410, v0401 - 47.3.2

          a[0]=1;
          b[0]=4;
          c[0]=0;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=0;
          b[2]=4;
          c[2]=1;

          a[3]=0;
          b[3]=4;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1400, v0311, v1301, v1310 - 47.3.4

          a[0]=1;
          b[0]=4;
          c[0]=0;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=1;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v1400, v0311, v1301, v0401 - 47.3.3

          a[0]=1;
          b[0]=4;
          c[0]=0;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=0;
          b[3]=4;
          c[3]=0;
        };
      };
    };        
  }
  else if (sC>(ValueType)3.0) //48
  {
    if(axis=='0') //48.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0140, v1031, v1130, v1040 - 48.1.1

          a[0]=0;
          b[0]=1;
          c[0]=4;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=4;
        }
        else
        {
          // v0140, v1031, v1130, v0131 - 48.1.3

          a[0]=0;
          b[0]=1;
          c[0]=4;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=3;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0140, v1031, v0041, v1040 - 48.1.4

          a[0]=0;
          b[0]=1;
          c[0]=4;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=0;
          b[2]=0;
          c[2]=4;

          a[3]=1;
          b[3]=0;
          c[3]=4;
        }
        else
        {
          // v0140, v1031, v0041, v0131 - 48.1.2

          a[0]=0;
          b[0]=1;
          c[0]=4;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=0;
          b[2]=0;
          c[2]=4;

          a[3]=0;
          b[3]=1;
          c[3]=3;
        };
      };
    }
    else if(axis=='1') //48.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1040, v0131, v1130, v0140 - 48.2.1

          a[0]=1;
          b[0]=0;
          c[0]=4;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=4;
        }
        else
        {
          // v1040, v0131, v1130, v1031 - 48.2.3

          a[0]=1;
          b[0]=0;
          c[0]=4;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1040, v0131, v0041, v0140 - 48.2.2

          a[0]=1;
          b[0]=0;
          c[0]=4;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=0;
          b[2]=0;
          c[2]=4;

          a[3]=0;
          b[3]=1;
          c[3]=4;
        }
        else
        {
          // v1040, v0131, v0041, v1031 - 48.2.4

          a[0]=1;
          b[0]=0;
          c[0]=4;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=0;
          b[2]=0;
          c[2]=4;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        };
      };
    }
    else //48.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1130, v0041, v0140, v1040 - 48.3.1

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=0;
          c[1]=4;

          a[2]=0;
          b[2]=1;
          c[2]=4;

          a[3]=1;
          b[3]=0;
          c[3]=4;
        }
        else
        {
          // v1130, v0041, v0140, v0131 - 48.3.2

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=0;
          c[1]=4;

          a[2]=0;
          b[2]=1;
          c[2]=4;

          a[3]=0;
          b[3]=1;
          c[3]=3;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1130, v0041, v1031, v1040 - 48.3.4

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=0;
          c[1]=4;

          a[2]=1;
          b[2]=0;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=4;
        }
        else
        {
          // v1130, v0041, v1031, v0131 - 48.3.3

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=0;
          c[1]=4;

          a[2]=1;
          b[2]=0;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=3;
        };
      };
    };        
  }
  else if (sD>(ValueType)3.0) //49
  {
    if(axis=='0') //49.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0113, v1004, v1103, v1013 - 49.1.1

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0113, v1004, v1103, v0104 - 49.1.3

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0113, v1004, v0014, v1013 - 49.1.4

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v0113, v1004, v0014, v0104 - 49.1.2

          a[0]=0;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //49.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1013, v0104, v1103, v0113 - 49.2.1

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1013, v0104, v1103, v1004 - 49.2.3

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1013, v0104, v0014, v0113 - 49.2.2

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1013, v0104, v0014, v1004 - 49.2.4

          a[0]=1;
          b[0]=0;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //49.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1103, v0014, v0113, v1013 - 49.3.1

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1103, v0014, v0113, v0104 - 49.3.2

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1103, v0014, v1004, v1013 - 49.3.4

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1103, v0014, v1004, v0104 - 49.3.3

          a[0]=1;
          b[0]=1;
          c[0]=0;

          a[1]=0;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=0;

          a[3]=0;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  }

  else if ((sA>(ValueType)2.0)&&(sB>(ValueType)2.0)) //5
  {
    // v3200, v2300, v2210, v2201 

    a[0]=3;
    b[0]=2;
    c[0]=0;

    a[1]=2;
    b[1]=3;
    c[1]=0;

    a[2]=2;
    b[2]=2;
    c[2]=1;

    a[3]=2;
    b[3]=2;
    c[3]=0;
  }
  else if ((sB>(ValueType)2.0)&&(sC>(ValueType)2.0)) //6
  {
    // v1220, v0320, v0230, v0221 

    a[0]=1;
    b[0]=2;
    c[0]=2;

    a[1]=0;
    b[1]=3;
    c[1]=2;

    a[2]=0;
    b[2]=2;
    c[2]=3;

    a[3]=0;
    b[3]=2;
    c[3]=2;
  }
  else if ((sC>(ValueType)2.0)&&(sD>(ValueType)2.0)) //7
  {
    // v1022, v0122, v0032, v0023 

    a[0]=1;
    b[0]=0;
    c[0]=2;

    a[1]=0;
    b[1]=1;
    c[1]=2;

    a[2]=0;
    b[2]=0;
    c[2]=3;

    a[3]=0;
    b[3]=0;
    c[3]=2;
  }
  else if ((sA>(ValueType)2.0)&&(sD>(ValueType)2.0)) //8
  {
    // v3002, v2102, v2012, v2003 

    a[0]=3;
    b[0]=0;
    c[0]=0;

    a[1]=2;
    b[1]=1;
    c[1]=0;

    a[2]=2;
    b[2]=0;
    c[2]=1;

    a[3]=2;
    b[3]=0;
    c[3]=0;
  }
  else if ((sA>(ValueType)2.0)&&(sC>(ValueType)2.0)) //9
  {
    // v3020, v2120, v2030, v2021 

    a[0]=3;
    b[0]=0;
    c[0]=2;

    a[1]=2;
    b[1]=1;
    c[1]=2;

    a[2]=2;
    b[2]=0;
    c[2]=3;

    a[3]=2;
    b[3]=0;
    c[3]=2;
  }
  else if ((sB>(ValueType)2.0)&&(sD>(ValueType)2.0)) //10
  {
    // v1202, v0302, v0212, v0203 

    a[0]=1;
    b[0]=2;
    c[0]=0;

    a[1]=0;
    b[1]=3;
    c[1]=0;

    a[2]=0;
    b[2]=2;
    c[2]=1;

    a[3]=0;
    b[3]=2;
    c[3]=0;
  }
  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)&&(sC>=(ValueType)1.0)) //23
  {
    // v3110, v2210, v2120, v2111 

    a[0]=3;
    b[0]=1;
    c[0]=1;

    a[1]=2;
    b[1]=2;
    c[1]=1;

    a[2]=2;
    b[2]=1;
    c[2]=2;

    a[3]=2;
    b[3]=1;
    c[3]=1;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //24
  {
    // v2210, v1310, v1220, v1211 

    a[0]=2;
    b[0]=2;
    c[0]=1;

    a[1]=1;
    b[1]=3;
    c[1]=1;

    a[2]=1;
    b[2]=2;
    c[2]=2;

    a[3]=1;
    b[3]=2;
    c[3]=1;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //25
  {
    // v2120, v1220, v1130, v1121 

    a[0]=2;
    b[0]=1;
    c[0]=2;

    a[1]=1;
    b[1]=2;
    c[1]=2;

    a[2]=1;
    b[2]=1;
    c[2]=3;

    a[3]=1;
    b[3]=1;
    c[3]=2;
  }

  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //26
  {
    // v1211, v0311, v0221, v0212 

    a[0]=1;
    b[0]=2;
    c[0]=1;

    a[1]=0;
    b[1]=3;
    c[1]=1;

    a[2]=0;
    b[2]=2;
    c[2]=2;

    a[3]=0;
    b[3]=2;
    c[3]=1;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //27
  {
    // v1121, v0221, v0131, v0122 

    a[0]=1;
    b[0]=1;
    c[0]=2;

    a[1]=0;
    b[1]=2;
    c[1]=2;

    a[2]=0;
    b[2]=1;
    c[2]=3;

    a[3]=0;
    b[3]=1;
    c[3]=2;
  }
  else if ((sB>=(ValueType)1.0)&&(sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //28
  {
    // v1112, v0212, v0122, v0113 

    a[0]=1;
    b[0]=1;
    c[0]=1;

    a[1]=0;
    b[1]=2;
    c[1]=1;

    a[2]=0;
    b[2]=1;
    c[2]=2;

    a[3]=0;
    b[3]=1;
    c[3]=1;
  }
  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //29
  {
    // v3101, v2201, v2111, v2102 

    a[0]=3;
    b[0]=1;
    c[0]=0;

    a[1]=2;
    b[1]=2;
    c[1]=0;

    a[2]=2;
    b[2]=1;
    c[2]=1;

    a[3]=2;
    b[3]=1;
    c[3]=0;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //30
  {
    // v2201, v1301, v1211, v1202 

    a[0]=2;
    b[0]=2;
    c[0]=0;

    a[1]=1;
    b[1]=3;
    c[1]=0;

    a[2]=1;
    b[2]=2;
    c[2]=1;

    a[3]=1;
    b[3]=2;
    c[3]=0;
  }
  else if ((sA>=(ValueType)1.0)&&(sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //31
  {
    // v2102, v1202, v1112, v1103 

    a[0]=2;
    b[0]=1;
    c[0]=0;

    a[1]=1;
    b[1]=2;
    c[1]=0;

    a[2]=1;
    b[2]=1;
    c[2]=1;

    a[3]=1;
    b[3]=1;
    c[3]=0;
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)&&(sD>=(ValueType)1.0)) //32
  {
    // v3011, v2111, v2021, v2012 

    a[0]=3;
    b[0]=0;
    c[0]=1;

    a[1]=2;
    b[1]=1;
    c[1]=1;

    a[2]=2;
    b[2]=0;
    c[2]=2;

    a[3]=2;
    b[3]=0;
    c[3]=1;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //33
  {
    // v2021, v1121, v1031, v1022 

    a[0]=2;
    b[0]=0;
    c[0]=2;

    a[1]=1;
    b[1]=1;
    c[1]=2;

    a[2]=1;
    b[2]=0;
    c[2]=3;

    a[3]=1;
    b[3]=0;
    c[3]=2;
  }
  else if ((sA>=(ValueType)1.0)&&(sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //34
  {
    // v2012, v1112, v1022, v1013 

    a[0]=2;
    b[0]=0;
    c[0]=1;

    a[1]=1;
    b[1]=1;
    c[1]=1;

    a[2]=1;
    b[2]=0;
    c[2]=2;

    a[3]=1;
    b[3]=0;
    c[3]=1;
  }

  else if ((sA>(ValueType)2.0)&&(sB>=(ValueType)1.0)) //50
  {
    if(axis=='0') //50.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2210, v3101, v3200, v3110 - 50.1.1

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=2;
          c[2]=0;

          a[3]=3;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2210, v3101, v3200, v2201   50.1.3

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2210, v3101, v2111, v3110   50.1.4

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=3;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2210, v3101, v2111, v2201   50.1.2

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=3;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //50.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3110, v2201, v3200, v2210   50.2.1

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=3;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v3110, v2201, v3200, v3101   50.2.3

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=3;
          b[2]=2;
          c[2]=0;

          a[3]=3;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3110, v2201, v2111, v2210   50.2.2

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v3110, v2201, v2111, v3101   50.2.4

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=3;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else //50.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3200, v2111, v2210, v3110   50.3.1

          a[0]=3;
          b[0]=2;
          c[0]=0;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=1;

          a[3]=3;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v3200, v2111, v2210, v2201   50.3.2

          a[0]=3;
          b[0]=2;
          c[0]=0;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=1;

          a[3]=2;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3200, v2111, v3101, v3110   50.3.4

          a[0]=3;
          b[0]=2;
          c[0]=0;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=3;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v3200, v2111, v3101, v2201   50.3.3

          a[0]=3;
          b[0]=2;
          c[0]=0;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=2;
          c[3]=0;
        };
      };
    };        
  }
  else if ((sB>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //51
  {
    if(axis=='0') //51.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0320, v1211, v1310, v1220   51.1.1

          a[0]=0;
          b[0]=3;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v0320, v1211, v1310, v0311   51.1.3

          a[0]=0;
          b[0]=3;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0320, v1211, v0221, v1220   51.1.4

          a[0]=0;
          b[0]=3;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=1;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v0320, v1211, v0221, v0311   51.1.2

          a[0]=0;
          b[0]=3;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=0;
          b[3]=3;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //51.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1220, v0311, v1310, v0320   51.2.1

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=2;
        }
        else
        {
          // v1220, v0311, v1310, v1211   51.2.3

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1220, v0311, v0221, v0320   51.2.2

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=0;
          b[3]=3;
          c[3]=2;
        }
        else
        {
          // v1220, v0311, v0221, v1211   51.2.4

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=3;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        };
      };
    }
    else //51.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1310, v0221, v0320, v1220   51.3.1

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=0;
          b[2]=3;
          c[2]=2;

          a[3]=1;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v1310, v0221, v0320, v0311   51.3.2

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=0;
          b[2]=3;
          c[2]=2;

          a[3]=0;
          b[3]=3;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1310, v0221, v1211, v1220   51.3.4

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v1310, v0221, v1211, v0311   51.3.3

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=1;
        };
      };
    };        
  }
  else if ((sC>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //52
  {
    if(axis=='0') //52.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0131, v1022, v1121, v1031   52.1.1

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v0131, v1022, v1121, v0122   52.1.3

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0131, v1022, v0032, v1031   52.1.4

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v0131, v1022, v0032, v0122   52.1.2

          a[0]=0;
          b[0]=1;
          c[0]=3;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      };
    }
    else if(axis=='1') //52.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1031, v0122, v1121, v0131   52.2.1

          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v1031, v0122, v1121, v1022   52.2.3

          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1031, v0122, v0032, v0131   52.2.2

          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v1031, v0122, v0032, v1022   52.2.4

          a[0]=1;
          b[0]=0;
          c[0]=3;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=0;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        };
      };
    }
    else //52.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1121, v0032, v0131, v1031   52.3.1

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v1121, v0032, v0131, v0122   52.3.2

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1121, v0032, v1022, v1031   52.3.4

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v1121, v0032, v1022, v0122   52.3.3

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //53
  {
    if(axis=='0') //53.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1112, v2003, v2102, v2012   53.1.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1112, v2003, v2102, v1103   53.1.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1112, v2003, v1013, v2012   53.1.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v1112, v2003, v1013, v1103   53.1.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //53.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2012, v1103, v2102, v1112   53.2.1

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2012, v1103, v2102, v2003   53.2.3

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2012, v1103, v1013, v1112   53.2.2

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2012, v1103, v1013, v2003   53.2.4

          a[0]=2;
          b[0]=0;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //53.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2102, v1013, v1112, v2012   53.3.1

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2102, v1013, v1112, v1103 - 53.3.2

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2102, v1013, v2003, v2012 - 53.3.4

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=0;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2102, v1013, v2003, v1103 - 53.3.3

          a[0]=2;
          b[0]=1;
          c[0]=0;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sB>(ValueType)2.0)) //54
  {
    if(axis=='0') //54.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1310, v2201, v2300, v2210 - 54.1.1

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=3;
          c[2]=0;

          a[3]=2;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1310, v2201, v2300, v1301 - 54.1.3

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=3;
          c[2]=0;

          a[3]=1;
          b[3]=3;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1310, v2201, v1211, v2210 - 54.1.4

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=2;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1310, v2201, v1211, v1301 - 54.1.2

          a[0]=1;
          b[0]=3;
          c[0]=1;

          a[1]=2;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=3;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //54.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2210, v1301, v2300, v1310 - 54.2.1

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=2;
          b[2]=3;
          c[2]=0;

          a[3]=1;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v2210, v1301, v2300, v2201 - 54.2.3

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=2;
          b[2]=3;
          c[2]=0;

          a[3]=2;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2210, v1301, v1211, v1310 - 54.2.2

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v2210, v1301, v1211, v2201 - 54.2.4

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=2;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else //54.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2300, v1211, v1310, v2210 - 54.3.1

          a[0]=2;
          b[0]=3;
          c[0]=0;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=1;

          a[3]=2;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v2300, v1211, v1310, v1301 - 54.3.2

          a[0]=2;
          b[0]=3;
          c[0]=0;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=3;
          c[2]=1;

          a[3]=1;
          b[3]=3;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2300, v1211, v2201, v2210 - 54.3.4

          a[0]=2;
          b[0]=3;
          c[0]=0;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v2300, v1211, v2201, v1301 - 54.3.3

          a[0]=2;
          b[0]=3;
          c[0]=0;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=3;
          c[3]=0;
        };
      };
    };        
  }
  else if ((sB>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //55
  {
    if(axis=='0') //55.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0230, v1121, v1220, v1130 - 55.1.1

          a[0]=0;
          b[0]=2;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v0230, v1121, v1220, v0221 - 55.1.3

          a[0]=0;
          b[0]=2;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0230, v1121, v0131, v1130 - 55.1.4

          a[0]=0;
          b[0]=2;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=1;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v0230, v1121, v0131, v0221 - 55.1.2

          a[0]=0;
          b[0]=2;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        };
      };
    }
    else if(axis=='1') //55.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1130, v0221, v1220, v0230    55.2.1

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=3;
        }
        else
        {
          // v1130, v0221, v1220, v1121   55.2.3

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1130, v0221, v0131, v0230   55.2.2

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=0;
          b[3]=2;
          c[3]=3;
        }
        else
        {
          // v1130, v0221, v0131, v1121   55.2.4

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=0;
          b[1]=2;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=3;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        };
      };
    }
    else //55.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1220, v0131, v0230, v1130   55.3.1

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=0;
          b[2]=2;
          c[2]=3;

          a[3]=1;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v1220, v0131, v0230, v0221   55.3.2

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=0;
          b[2]=2;
          c[2]=3;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1220, v0131, v1121, v1130   55.3.4

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v1220, v0131, v1121, v0221   55.3.3

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        };
      };
    };        
  }
  else if ((sC>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //56
  {
    if(axis=='0') //56.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0122, v1013, v1112, v1022   56.1.1

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v0122, v1013, v1112, v0113   56.1.3

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0122, v1013, v0023, v1022   56.1.4

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v0122, v1013, v0023, v0113   56.1.2

          a[0]=0;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //56.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)1.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1022, v0113, v1112, v0122   56.2.1

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1022, v0113, v1112, v1013   56.2.3

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1022, v0113, v0023, v0122   56.2.2

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1022, v0113, v0023, v1013   56.2.4

          a[0]=1;
          b[0]=0;
          c[0]=2;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=1;
        };
      };
    }
    else //56.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1112, v0023, v0122, v1022   56.3.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1112, v0023, v0122, v0113   56.3.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1112, v0023, v1013, v1022   56.3.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1112, v0023, v1013, v0113   56.3.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=1;

          a[3]=0;
          b[3]=1;
          c[3]=1;
        };
      };
    };        
  }
  else if ((sA>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //57
  {
    if(axis=='0') //57.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2111, v3002, v3101, v3011   57.1.1

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2111, v3002, v3101, v2102   57.1.3

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2111, v3002, v2012, v3011   57.1.4

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v2111, v3002, v2012, v2102   57.1.2

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=3;
          b[1]=0;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //57.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3011, v2102, v3101, v2111   57.2.1

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v3011, v2102, v3101, v3002   57.2.3

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=3;
          b[2]=1;
          c[2]=0;

          a[3]=3;
          b[3]=0;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3011, v2102, v2012, v2111   57.2.2

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v3011, v2102, v2012, v3002   57.2.4

          a[0]=3;
          b[0]=0;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=0;
        };
      };
    }
    else //57.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)1.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3101, v2012, v2111, v3011   57.3.1

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v3101, v2012, v2111, v2102   57.3.2

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3101, v2012, v3002, v3011   57.3.4

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=0;
          c[2]=0;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        }
        else
        {
          // v3101, v2012, v3002, v2102   57.3.3

          a[0]=3;
          b[0]=1;
          c[0]=0;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=0;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      };
    };        
  }
  else if ((sA>(ValueType)2.0)&&(sC>=(ValueType)1.0)) //58
  {
    if(axis=='0') //58.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2120, v3011, v3110, v3020   58.1.1

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v2120, v3011, v3110, v2111   58.1.3

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2120, v3011, v2021, v3020   58.1.4

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=2;

          a[3]=3;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v2120, v3011, v2021, v2111   58.1.2

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=3;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=2;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //58.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3020, v2111, v3110, v2120   58.2.1

          a[0]=3;
          b[0]=0;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v3020, v2111, v3110, v3011   58.2.3

          a[0]=3;
          b[0]=0;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=3;
          b[2]=1;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3020, v2111, v2021, v2120   58.2.2

          a[0]=3;
          b[0]=0;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=2;

          a[3]=2;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v3020, v2111, v2021, v3011   58.2.4

          a[0]=3;
          b[0]=0;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=0;
          c[2]=2;

          a[3]=3;
          b[3]=0;
          c[3]=1;
        };
      };
    }
    else //58.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v3110, v2021, v2120, v3020   58.3.1

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=2;

          a[3]=3;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v3110, v2021, v2120, v2111   58.3.2

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=2;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v3110, v2021, v3011, v3020   58.3.4

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=3;
          b[2]=0;
          c[2]=1;

          a[3]=3;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v3110, v2021, v3011, v2111   58.3.3

          a[0]=3;
          b[0]=1;
          c[0]=1;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=3;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        };
      };
    };        
  }
  else if ((sB>(ValueType)2.0)&&(sD>=(ValueType)1.0)) //59
  {
    if(axis=='0') //59.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0311, v1202, v1301, v1211   59.1.1

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v0311, v1202, v1301, v0302   59.1.3

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0311, v1202, v0212, v1211   59.1.4

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v0311, v1202, v0212, v0302   59.1.2

          a[0]=0;
          b[0]=3;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //59.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1211, v0302, v1301, v0311   59.2.1

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=0;
          b[3]=3;
          c[3]=1;
        }
        else
        {
          // v1211, v0302, v1301, v1202   59.2.3

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=1;
          b[2]=3;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1211, v0302, v0212, v0311   59.2.2

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=1;
       }
        else
        {
          // v1211, v0302, v0212, v1202   59.2.4

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=3;
          c[1]=0;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else //59.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1301, v0212, v0311, v1211   59.3.1

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1301, v0212, v0311, v0302   59.3.2

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=3;
          c[2]=1;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1301, v0212, v1202, v1211   59.3.4

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1301, v0212, v1202, v0302   59.3.3

          a[0]=1;
          b[0]=3;
          c[0]=0;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=3;
          c[3]=0;
        };
      };
    };        
  }
  else if ((sA>=(ValueType)1.0)&&(sC>(ValueType)2.0)) //60
  {
    if(axis=='0') //60.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1130, v2021, v2120, v2030   60.1.1

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v1130, v2021, v2120, v1121   60.1.3

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1130, v2021, v1031, v2030   60.1.4

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=3;

          a[3]=2;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v1130, v2021, v1031, v1121   60.1.2

          a[0]=1;
          b[0]=1;
          c[0]=3;

          a[1]=2;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=3;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        };
      };
    }
    else if(axis=='1') //60.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2030, v1121, v2120, v1130   60.2.1

          a[0]=2;
          b[0]=0;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v2030, v1121, v2120, v2021   60.2.3

          a[0]=2;
          b[0]=0;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2030, v1121, v1031, v1130   60.2.2

          a[0]=2;
          b[0]=0;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=3;

          a[3]=1;
          b[3]=1;
          c[3]=3;
        }
        else
        {
          // v2030, v1121, v1031, v2021   60.2.4

          a[0]=2;
          b[0]=0;
          c[0]=3;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=0;
          c[2]=3;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        };
      };
    }
    else //60.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2120, v1031, v1130, v2030   60.3.1

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=3;

          a[3]=2;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v2120, v1031, v1130, v1121   60.3.2

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=1;
          b[2]=1;
          c[2]=3;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2120, v1031, v2021, v2030   60.3.4

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=2;
          b[2]=0;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=3;
        }
        else
        {
          // v2120, v1031, v2021, v1121   60.3.3

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=0;
          c[1]=3;

          a[2]=2;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        };
      };
    };        
  }
  else if ((sB>=(ValueType)1.0)&&(sD>(ValueType)2.0)) //61
  {
    if(axis=='0') //61.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0212, v1103, v1202, v1112   61.1.1

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v0212, v1103, v1202, v0203   61.1.3

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0212, v1103, v0113, v1112   61.1.4

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v0212, v1103, v0113, v0203   61.1.2

          a[0]=0;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //61.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1112, v0203, v1202, v0212   61.2.1

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1112, v0203, v1202, v1103   61.2.3

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1112, v0203, v0113, v0212   61.2.2

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v1112, v0203, v0113, v1103   61.2.4

          a[0]=1;
          b[0]=1;
          c[0]=1;

          a[1]=0;
          b[1]=2;
          c[1]=0;

          a[2]=0;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else //61.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)1.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1202, v0113, v0212, v1112   61.3.1

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1202, v0113, v0212, v0203   61.3.2

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1202, v0113, v1103, v1112   61.3.4

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1202, v0113, v1103, v0203   61.3.3

          a[0]=1;
          b[0]=2;
          c[0]=0;

          a[1]=0;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=0;

          a[3]=0;
          b[3]=2;
          c[3]=0;
        };
      };
    };        
  }

  else if (sA>(ValueType)2.0) //36
  {
    // v2111, v3011, v3101, v3110 

    a[0]=2;
    b[0]=1;
    c[0]=1;

    a[1]=3;
    b[1]=0;
    c[1]=1;

    a[2]=3;
    b[2]=1;
    c[2]=0;

    a[3]=3;
    b[3]=1;
    c[3]=1;
  }
  else if (sB>(ValueType)2.0) //37
  {
    // v0311, v1211, v1301, v1310 

    a[0]=0;
    b[0]=3;
    c[0]=1;

    a[1]=1;
    b[1]=2;
    c[1]=1;

    a[2]=1;
    b[2]=3;
    c[2]=0;

    a[3]=1;
    b[3]=3;
    c[3]=1;
  }
  else if (sC>(ValueType)2.0) //38
  {
    // v0131, v1031, v1121, v1130 

    a[0]=0;
    b[0]=1;
    c[0]=3;

    a[1]=1;
    b[1]=0;
    c[1]=3;

    a[2]=1;
    b[2]=1;
    c[2]=2;

    a[3]=1;
    b[3]=1;
    c[3]=3;
  }
  else if (sD>(ValueType)2.0) //39
  {
    // v0113, v1013, v1103, v1112 

    a[0]=0;
    b[0]=1;
    c[0]=1;

    a[1]=1;
    b[1]=0;
    c[1]=1;

    a[2]=1;
    b[2]=1;
    c[2]=0;

    a[3]=1;
    b[3]=1;
    c[3]=1;
  }

  else if ((sC<(ValueType)1.0)&&(sD<(ValueType)1.0)) //40
  {
    // v1211, v2111, v2201, v2210 

    a[0]=1;
    b[0]=2;
    c[0]=1;

    a[1]=2;
    b[1]=1;
    c[1]=1;

    a[2]=2;
    b[2]=2;
    c[2]=0;

    a[3]=2;
    b[3]=2;
    c[3]=1;
  }
  else if ((sA<(ValueType)1.0)&&(sD<(ValueType)1.0)) //41
  {
    // v0221, v1121, v1211, v1220 

    a[0]=0;
    b[0]=2;
    c[0]=2;

    a[1]=1;
    b[1]=1;
    c[1]=2;

    a[2]=1;
    b[2]=2;
    c[2]=1;

    a[3]=1;
    b[3]=2;
    c[3]=2;
  }
  else if ((sA<(ValueType)1.0)&&(sB<(ValueType)1.0)) //42
  {
    // v0122, v1022, v1112, v1121 

    a[0]=0;
    b[0]=1;
    c[0]=2;

    a[1]=1;
    b[1]=0;
    c[1]=2;

    a[2]=1;
    b[2]=1;
    c[2]=1;

    a[3]=1;
    b[3]=1;
    c[3]=2;
  }
  else if ((sB<(ValueType)1.0)&&(sC<(ValueType)1.0)) //43
  {
    // v1112, v2012, v2102, v2111 

    a[0]=1;
    b[0]=1;
    c[0]=1;

    a[1]=2;
    b[1]=0;
    c[1]=1;

    a[2]=2;
    b[2]=1;
    c[2]=0;

    a[3]=2;
    b[3]=1;
    c[3]=1;
  }
  else if ((sB<(ValueType)1.0)&&(sD<(ValueType)1.0)) //44
  {
    // v1121, v2021, v2111, v2120 

    a[0]=1;
    b[0]=1;
    c[0]=2;

    a[1]=2;
    b[1]=0;
    c[1]=2;

    a[2]=2;
    b[2]=1;
    c[2]=1;

    a[3]=2;
    b[3]=1;
    c[3]=2;
  }
  else if ((sA<(ValueType)1.0)&&(sC<(ValueType)1.0)) //45
  {
    // v0212, v1112, v1202, v1211 

    a[0]=0;
    b[0]=2;
    c[0]=1;

    a[1]=1;
    b[1]=1;
    c[1]=1;

    a[2]=1;
    b[2]=2;
    c[2]=0;

    a[3]=1;
    b[3]=2;
    c[3]=1;
  }

  else if (sA<(ValueType)1.0) //62
  {
    if(axis=='0') //62.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v0221, v1112, v1211, v1121   62.1.1

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v0221, v1112, v1211, v0212   62.1.3

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v0221, v1112, v0122, v1121   62.1.4

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v0221, v1112, v0122, v0212   62.1.2

          a[0]=0;
          b[0]=2;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //62.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1121, v0212, v1211, v0221   62.2.1

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v1121, v0212, v1211, v1112   62.2.3

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1121, v0212, v0122, v0221   62.2.2

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v1121, v0212, v0122, v1112   62.2.4

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=0;
          b[1]=2;
          c[1]=1;

          a[2]=0;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else //62.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1211, v0122, v0221, v1121   62.3.1

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1211, v0122, v0221, v0212   62.3.2

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=0;
          b[2]=2;
          c[2]=2;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1211, v0122, v1112, v1121   62.3.4

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1211, v0122, v1112, v0212   62.3.3

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=0;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=0;
          b[3]=2;
          c[3]=1;
        };
      };
    };        
  }
  else if (sB<(ValueType)1.0) //63
  {
    if(axis=='0') //63.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1121, v2012, v2111, v2021   63.1.1

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1121, v2012, v2111, v1112   63.1.3

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1121, v2012, v1022, v2021   63.1.4

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v1121, v2012, v1022, v1112   63.1.2

          a[0]=1;
          b[0]=1;
          c[0]=2;

          a[1]=2;
          b[1]=0;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //63.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0+(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2021, v1112, v2111, v1121   63.2.1

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v2021, v1112, v2111, v2012   63.2.3

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2021, v1112, v1022, v1121   63.2.2

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v2021, v1112, v1022, v2012   63.2.4

          a[0]=2;
          b[0]=0;
          c[0]=2;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=0;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=1;
        };
      };
    }
    else //63.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2111, v1022, v1121, v2021   63.3.1

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v2111, v1022, v1121, v1112   63.3.2

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2111, v1022, v2012, v2021   63.3.4

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=2;
          b[3]=0;
          c[3]=2;
        }
        else
        {
          // v2111, v1022, v2012, v1112   63.3.3

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=0;
          c[1]=2;

          a[2]=2;
          b[2]=0;
          c[2]=1;

          a[3]=1;
          b[3]=1;
          c[3]=1;
        };
      };
    };        
  }
  else if (sC<(ValueType)1.0) //64
  {
    if(axis=='0') //64.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1211, v2102, v2201, v2111   64.1.1

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1211, v2102, v2201, v1202   64.1.3

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1211, v2102, v1112, v2111   64.1.4

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v1211, v2102, v1112, v1202   64.1.2

          a[0]=1;
          b[0]=2;
          c[0]=1;

          a[1]=2;
          b[1]=1;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      };
    }
    else if(axis=='1') //64.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2111, v1202, v2201, v1211   64.2.1

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v2111, v1202, v2201, v2102  64.2.3

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=2;
          b[2]=2;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2111, v1202, v1112, v1211   64.2.2

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        }
        else
        {
          // v2111, v1202, v1112, v2102   64.2.4

          a[0]=2;
          b[0]=1;
          c[0]=1;

          a[1]=1;
          b[1]=2;
          c[1]=0;

          a[2]=1;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=0;
        };
      };
    }
    else //64.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0+(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0+(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2201, v1112, v1211, v2111   64.3.1

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2201, v1112, v1211, v1202   64.3.2

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2201, v1112, v2102, v2111   64.3.4

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        }
        else
        {
          // v2201, v1112, v2102, v1202   64.3.3

          a[0]=2;
          b[0]=2;
          c[0]=0;

          a[1]=1;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=1;
          c[2]=0;

          a[3]=1;
          b[3]=2;
          c[3]=0;
        };
      };
    };        
  }
  else if (sD<(ValueType)1.0) //65
  {
    if(axis=='0') //65.1
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v1220, v2111, v2210, v2120   65.1.1

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1220, v2111, v2210, v1211   65.1.3

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v1220, v2111, v1121, v2120   65.1.4

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=2;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v1220, v2111, v1121, v1211    65.1.2

          a[0]=1;
          b[0]=2;
          c[0]=2;

          a[1]=2;
          b[1]=1;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        };
      };
    }
    else if(axis=='1') //65.2
    {
      s1=(sA+sB-sC-sD)/(ValueType)2.0-(ValueType)0.5;
      s2=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2120, v1211, v2210, v1220   65.2.1

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v2120, v1211, v2210, v2111    65.2.3

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=2;
          b[2]=2;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2120, v1211, v1121, v1220   65.2.2

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=1;
          b[3]=2;
          c[3]=2;
        }
        else
        {
          // v2120, v1211, v1121, v2111   65.2.4

          a[0]=2;
          b[0]=1;
          c[0]=2;

          a[1]=1;
          b[1]=2;
          c[1]=1;

          a[2]=1;
          b[2]=1;
          c[2]=2;

          a[3]=2;
          b[3]=1;
          c[3]=1;
        };
      };
    }
    else //65.3
    {
      s1=(sB+sC-sD-sA)/(ValueType)2.0-(ValueType)0.5;
      s2=(sA+sC-sB-sD)/(ValueType)2.0-(ValueType)0.5;
      if (s1>(ValueType)0)
      {
        if (s2>(ValueType)0)
        {
          // v2210, v1121, v1220, v2120   65.3.1

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=2;

          a[3]=2;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v2210, v1121, v1220, v1211   65.3.2

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=1;
          b[2]=2;
          c[2]=2;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        };
      }
      else
      {
        if (s2>(ValueType)0)
        {
          // v2210, v1121, v2111, v2120   65.3.4

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=2;
          b[3]=1;
          c[3]=2;
        }
        else
        {
          // v2210, v1121, v2111, v1211   65.3.3

          a[0]=2;
          b[0]=2;
          c[0]=1;

          a[1]=1;
          b[1]=1;
          c[1]=2;

          a[2]=2;
          b[2]=1;
          c[2]=1;

          a[3]=1;
          b[3]=2;
          c[3]=1;
        };
      };
    };        
  }

  else //35
  {
    // v2111, v1211, v1121, v1112 

    a[0]=2;
    b[0]=1;
    c[0]=1;

    a[1]=1;
    b[1]=2;
    c[1]=1;

    a[2]=1;
    b[2]=1;
    c[2]=2;

    a[3]=1;
    b[3]=1;
    c[3]=1;
  };

  for(i=0; (i<4); i++)
  {
    Coord[i*3  ]=xD+(xA-xD)*a[i]+(xB-xD)*b[i]+(xC-xD)*c[i];
    Coord[i*3+1]=yD+(yA-yD)*a[i]+(yB-yD)*b[i]+(yC-yD)*c[i];
    Coord[i*3+2]=zD+(zA-zD)*a[i]+(zB-zD)*b[i]+(zC-zD)*c[i];
  };

  return;
}


#endif

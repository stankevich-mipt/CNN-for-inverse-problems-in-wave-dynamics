#ifndef InterpolationHeaderAleanera2012
#define InterpolationHeaderAleanera2012

#include "math.h"

///////////////////////////FOR SWIFT///////////////////////////////////////////////////

////////////////////////Interpolation/////////////////////////////////////////

template <typename ValueType>
ValueType Interpolate_3D(
//номера вершин тетраэдра ABCD во всей сетке: a<b<c<d

const int N, 
//порядок интерполяции

const char How, 
//тип интерполяции:
//'L' - линейная, 
//'M' - монотонная, 
//'P' - полиномом

const char Axis,
//Номер оси для кусочно-линейной интерполяции и интерполяции с ограничителем. '0', '1', '2', остальное равносильно '2'

const ValueType * R,
//массив из 15 координат
//x, y, z векторов R, Ra, Rb, Rc, Rd
//x, y, z, xA, yA, zA, xB, yB, zB, xC, ...

const ValueType * V);
//массив из значений в (N+1)*(N+2)*(N+3)/6 опорных точках
//точки вершин: A, B, C, D - по 1 значению в каждой
//точки ребер: AB, AC, AD, BC, BD, CD - по N-1 значений для каждой
//точки граней: ABC, ABD, ACD, BCD - по (N-1)*(N-2)/2 значений для каждой
//точки объема: ABCD - (N-1)*(N-2)*(N-3)/6 значений
//вначале даются значения в вершинах, зачем значения во всеx точках стороны AB, 
//затем значения во всех точках стороны AC и т.д.

////////////////Координаты вершин маленького тетраэдра////////////////

template <typename ValueType>
void Nodes_Coord_3D(
ValueType * Coord,
const int N, 
const char Axis,
const ValueType * R);
//Возвращает координаты вершин малого тетраэдра, в который попала точка, Coord - массив из 12, как обычно x1, y1, z1, x2, y2, z2, x3, ...

template <typename ValueType>
void Nodes_Coord_3D_2(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD);

template <typename ValueType>
void Nodes_Coord_3D_3(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD);

template <typename ValueType>
void Nodes_Coord_3D_4(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD);

/* void Nodes_Coord_3D_5(ValueType * Coord, const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD); */

///////////////////Перестановки для контактов//////////////////

int GetRelativeOrientationSide(int * firstNodes, int * secondNodes);
// firstNodes - 3 номера вершин 1го треугольника, 
//secondNodes - 3 номера вершин 2го треугольника

int GetCorrespondingIndexSide(int relativeOrientation, int firstIndex, int number);
// relativeOrientation число до 1210, характеризующее относительную ориентацию; 
//firstIndex индекс в 1ом треугольнике, которому нужно найти соответствие во 2ом;
//number - порядок интерполяции

int GetRelativeOrientationEdge(int * firstNodes, int * secondNodes);
// firstNodes - 2 номера вершин 1го отрезка, 
//secondNodes - 2 номера вершин 2го отрезка

int GetCorrespondingIndexEdge(int relativeOrientation, int firstIndex, int number);
// relativeOrientation число (0 или 1), характеризующее относительную ориентацию; 
//firstIndex индекс в 1ом отрезке, которому нужно найти соответствие во 2ом;
//number - порядок интерполяции

//////////////функции для нахождения радиусов опорных точек//////////////

template <typename ValueType>
void Interpolation_3D_Get_Edge(
//возвращает координату опорной точки из массива на ребре AB, a<b
ValueType * Result, 
//массив из 3 координат опорной точки, которые надо найти
const int N, 
//порядок интерполяции
const ValueType * R, 
//массив из 6 координат: xA, yA, zA, xB, yB, zB
const int I);
//номер опорной точки в массиве на данном ребрe A: 0<=I<N-1

template <typename ValueType>
void Interpolation_3D_Get_Side(
//возвращает координату опорной точки из массива на грани ABC, a<b<c
ValueType * Result, 
//массив из 3 координат опорной точки, которые надо найти
const int N, 
//порядок интерполяции
const ValueType * R, 
//массив из 9 координат: xA, yA, zA, xB, yB, zB, xC, yC, zC
const int I);
//номер опорной точки в массиве на данной грани ABC: 0<=I<(N-1)*(N-2)/2

template <typename ValueType>
void Interpolation_3D_Get_Volume(
//возвращает координату опорной точки из массива внутри тетраэдра ABCD, a<b<c<d
ValueType * Result, 
//массив из 3 координат опорной точки, которые надо найти
const int N, 
//порядок интерполяции
const ValueType * R, 
//массив из 12 координат: xA, yA, zA, xB, yB, zB, xC, yC, zC, xD, yD, zD
const int I);
//номер опорной точки в массиве внутри данного тетраэдра ABCD: 0<=I<(N-1)*(N-2)*(N-3)/6

///////////////////////Permitaions////////////////////////////////////////////

template <typename ValueType>
void InterpolationPermitation_Edge(const int N, const int * P, ValueType * V);
//N>1, P[2], V[N-1]

template <typename ValueType>
void InterpolationPermitation_Side(const int N, const int * P, ValueType * V);
//N>2, P[3], V[(N-1)*(N-2)/2]

template <typename ValueType>
void InterpolationPermitation_Volume(const int N, const int * P, ValueType * V);
//N>3, P[4], V[(N-1)*(N-2)*(N-3)/6]

///////////////////////Courant////////////////////////////////////////////////

template <typename ValueType>
ValueType Interpolation_MinAltitude_3D_Swift(char * Axis, const int N, const ValueType * R);
//R -- массив из 12 чисел с координатами вершин тетраэдра. xA, yA, zA, xB, ...

///////////////////////MaxAltitude////////////////////////////////

template <typename ValueType>
ValueType Interpolation_MaxAltitude_3D_Swift(char Axis, const int N, const ValueType * R);
//R -- массив из 12 чисел с координатами вершин тетраэдра. xA, yA, zA, xB, ...

/////////////////////////////////interpolation_2_4///////////////////////////////////////////////////////////

template <typename ValueType>
ValueType Interpolate_3D_interpolation_2_4(
const char How, 
const char Axis,
const ValueType * R,
const ValueType * V);

template <typename ValueType>
ValueType Interpolation_MinAltitude_3D_Swift_interpolation_2_4(char * Axis, const ValueType * R);



///////////////////////OTHER///////////////////////////////////////////////////////////

/////////////////////////////////interpolation_2_4///////////////////////////////////////////////////////////

template <typename ValueType>
ValueType interpolation_MinAltitude_3D_interpolation_2_4_axis_1(
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD);

template <typename ValueType>
void interpolation_ABCD_Coordinate(const int a, const int b, const int c, const int d, const int N,
                                   ValueType * X, ValueType * Y, ValueType * Z, 
                                   const ValueType xA, const ValueType yA, const ValueType zA,  
                                   const ValueType xB, const ValueType yB, const ValueType zB,  
                                   const ValueType xC, const ValueType yC, const ValueType zC,  
                                   const ValueType xD, const ValueType yD, const ValueType zD);


template <typename ValueType>
ValueType interpolationAreaOfTriangle_3D(
ValueType xA, ValueType yA, ValueType zA,
ValueType xB, ValueType yB, ValueType zB,
ValueType xC, ValueType yC, ValueType zC);

template <typename ValueType>
ValueType interpolation_Max_2(const ValueType s1, const ValueType s2);

template <typename ValueType>
ValueType interpolation_MinAndNumber_3(char * a,
const ValueType s0, const ValueType s1, const ValueType s2);

template <typename ValueType>
ValueType interpolation_Max_4(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4);

template <typename ValueType>
ValueType interpolation_Max_5(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4, const ValueType s5);

template <typename ValueType>
ValueType interpolation_MinAltitude_3D(char * Axis, const int N,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD);

template <typename ValueType>
ValueType interpolation_Min_4(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4);

template <typename ValueType>
ValueType interpolation_Min_6(const ValueType s1, const ValueType s2, 
const ValueType s3, const ValueType s4, const ValueType s5, ValueType s6);

//////////////////////////////////////////////N=1/////////////////////////////////////////

//2D

template <typename ValueType>
ValueType Interpolate_2D_Linear(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v100, const ValueType v010, const ValueType v001);

//3D

template <typename ValueType>
ValueType Interpolate_3D_Linear(
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v1000, const ValueType v0100, const ValueType v0010, const ValueType v0001);

//////////////////////////////////////////////N=2//////////////////////////////////////////

//2D

template <typename ValueType>
ValueType Interpolate_2D_Power_2(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v200, const ValueType v020, const ValueType v002,
const ValueType v110, const ValueType v011, const ValueType v101);

template <typename ValueType>
ValueType Interpolate_2D_Linear_2(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v200, const ValueType v020, const ValueType v002,
const ValueType v110, const ValueType v011, const ValueType v101);

template <typename ValueType>
ValueType Interpolate_2D_Mono_2(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v200, const ValueType v020, const ValueType v002,
const ValueType v110, const ValueType v011, const ValueType v101);

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
const ValueType v1001, const ValueType v1010, const ValueType v0101);

template <typename ValueType>
ValueType Interpolate_3D_Linear_2(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v2000, const ValueType v0200, const ValueType v0020, const ValueType v0002,
const ValueType v1100, const ValueType v0110, const ValueType v0011, 
const ValueType v1001, const ValueType v1010, const ValueType v0101);

template <typename ValueType>
ValueType Interpolate_3D_Mono_2(const char axis,
const ValueType x, const ValueType y, const ValueType z,
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC,
const ValueType xD, const ValueType yD, const ValueType zD,
const ValueType v2000, const ValueType v0200, const ValueType v0020, const ValueType v0002,
const ValueType v1100, const ValueType v0110, const ValueType v0011, 
const ValueType v1001, const ValueType v1010, const ValueType v0101);

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
const ValueType v111);

template <typename ValueType>
ValueType Interpolate_2D_Linear_3(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v300, const ValueType v030, const ValueType v003,
const ValueType v210, const ValueType v021, const ValueType v102, 
const ValueType v120, const ValueType v012, const ValueType v201,
const ValueType v111);

template <typename ValueType>
ValueType Interpolate_2D_Mono_3(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
const ValueType v300, const ValueType v030, const ValueType v003,
const ValueType v210, const ValueType v021, const ValueType v102,
const ValueType v120, const ValueType v012, const ValueType v201,
const ValueType v111);

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
const ValueType v1110, const ValueType v0111, const ValueType v1011, const ValueType v1101);

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
const ValueType v1110, const ValueType v0111, const ValueType v1011, const ValueType v1101);

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
const ValueType v1110, const ValueType v0111, const ValueType v1011, const ValueType v1101);

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
const ValueType v211, const ValueType v121, const ValueType v112);

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
const ValueType v211, const ValueType v121, const ValueType v112);

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
const ValueType v211, const ValueType v121, const ValueType v112);

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
const ValueType v1111);

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
const ValueType v1111);

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
const ValueType v1111);

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
const ValueType v221, const ValueType v122, const ValueType v212);

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
const ValueType v221, const ValueType v122, const ValueType v212);

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
const ValueType v221, const ValueType v122, const ValueType v212);

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
const ValueType v2111, const ValueType v1211, const ValueType v1121, const ValueType v1112);

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
const ValueType v2111, const ValueType v1211, const ValueType v1121, const ValueType v1112);

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
const ValueType v2111, const ValueType v1211, const ValueType v1121, const ValueType v1112);

//////////////////
template <typename ValueType>
ValueType interpolationVectorProduct2D(
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB);

template <typename ValueType>
ValueType interpolationCompositionalProduct3D(
const ValueType xA, const ValueType yA, const ValueType zA,
const ValueType xB, const ValueType yB, const ValueType zB,
const ValueType xC, const ValueType yC, const ValueType zC);

template <typename ValueType>
void interpolationInitial_2D(
const ValueType x, const ValueType y,
const ValueType xA, const ValueType yA,
const ValueType xB, const ValueType yB,
const ValueType xC, const ValueType yC,
ValueType * sA,
ValueType * sB,
ValueType * sC);

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
ValueType * sD);

template <typename ValueType>
void interpolationMinMax_2D(ValueType * Min, ValueType * Max, 
const ValueType A, const ValueType B, const ValueType C);

template <typename ValueType>
void interpolationMinMax_3D(ValueType * Min, ValueType * Max, 
const ValueType A, const ValueType B, const ValueType C, const ValueType D);

#endif

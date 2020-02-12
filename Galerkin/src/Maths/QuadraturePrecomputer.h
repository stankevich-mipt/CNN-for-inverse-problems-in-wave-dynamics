#pragma once

#include "../../../3rdparty/quadrature_integration/legendre_rule.h"
#include "../../../3rdparty/quadrature_integration/triangle_fekete_rule.hpp"
#include "../../../3rdparty/quadrature_integration/tetrahedron_arbq_rule.hpp"

#include "Spaces.h"
#include <vector>
#include <assert.h>

struct QuadraturePrecomputer
{
  template <typename Space>
  static void BuildQuadrature(int order, std::vector<typename Space::Scalar>& weights, std::vector<typename Space::Vector>& points);
};

template <>
void QuadraturePrecomputer::BuildQuadrature<Space1>(int order, std::vector<Space1::Scalar>& weights, std::vector<Space1::Vector>& points)
{
  typedef Space1::Scalar Scalar;
  // http://people.sc.fsu.edu/~jburkardt/cpp_src/legendre_rule/legendre_rule.html
  // if we have N Gauss points then integral is absoltely presice for polynoms of up to 2N-1 degree

  weights.resize(order);
  points.resize(order); 
  const int    kind = 1; // gauss
  const Scalar alpha = 0;
  const Scalar beta = 0;
  const int numberOfGaussPoints = int((order + Scalar(1.0)) / 2 + 1 - std::numeric_limits<Scalar>::epsilon());
  cgqf(numberOfGaussPoints , kind, alpha, beta, 0, 1, points.data(), weights.data());
}

template <>
void QuadraturePrecomputer::BuildQuadrature<Space2>(int order, std::vector<Space2::Scalar>& weights, std::vector<Space2::Vector>& points)
{
  typedef Space2::Scalar Scalar;

  // variable "rule" sets number of gauss-legendre-lobatto points, i.e. accuracy of integration
  // http://people.sc.fsu.edu/~jburkardt/cpp_src/triangle_fekete_rule/triangle_fekete_rule.html

  /*Rule	Order	Precision
     1	   10	        3
     2	   28	        6
     3	   55	        9
     4	   91	       12
     5	   91	       12
     6	  136	       15
     7	  190	       18 */

  const int Precisions[] = {3, 6, 9, 12, 12, 15, 18};
  int rule = 1;
  while (rule <= 7 && order > Precisions[rule - 1])
    rule++;

  int rule_num = ::fekete_rule_num();
  assert(rule <= rule_num);

  int order_num = ::fekete_order_num(rule);

  std::vector<Scalar> xy(2 * order_num);
  std::vector<Scalar> w(order_num);

  ::fekete_rule(rule, order_num, xy.data(), w.data());

  for (int i = 0; i < order_num; ++i)
  {
    weights.push_back(w[i] * Scalar(0.5) /* area of unit triangle */);
    points.push_back(Space2::Vector(xy[2 * i], xy[2 * i + 1]));
  }
}

template <>
void QuadraturePrecomputer::BuildQuadrature<Space3>(int order, std::vector<Space3::Scalar>& weights, std::vector<Space3::Vector>& points)
{
  typedef Space3::Scalar Scalar;

  // from http://people.sc.fsu.edu/~jburkardt/c_src/tetrahedron_arbq_rule/tetrahedron_arbq_rule.html

  int order_num = ::tetrahedron_arbq_size(order);

  std::vector<Scalar> xyz(3 * order_num);
  std::vector<Scalar> w(order_num);

  ::tetrahedron_arbq(order, order_num, xyz.data(), w.data());

  for (int i = 0; i < order_num; ++i)
  {
    weights.push_back(w[i] / (4 * sqrt(2.0))); 
    /* 
       volume of their reference tetrahedron = sqrt(8)/3  
       and 
       volume of our unit tetrahedron = 1/6
    */
    Scalar* uvw = ::ref_to_koorn(xyz.data() + 3 * i);
    points.push_back(Space3::Vector(uvw[0] + Scalar(1.0), uvw[1] + Scalar(1.0), uvw[2] + Scalar(1.0)) * Scalar(0.5));
    free(uvw);
  }
}
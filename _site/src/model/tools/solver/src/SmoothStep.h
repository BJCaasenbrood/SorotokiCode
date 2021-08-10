#ifndef SMOOTHSTEP_H
#define SMOOTHSTEP_H

#include <algorithm>
#include <cmath>
#include "Eigen/Dense"
using namespace Eigen;

using namespace std;
typedef Eigen::VectorXd Vxd;
typedef Eigen::MatrixXd Mxd;

float clamp(double x, double lowerlimit, double upperlimit) {
  if (x < lowerlimit)
    x = lowerlimit;
  if (x > upperlimit)
    x = upperlimit;
  return x;
}

Vxd clampmag(Vxd x, double d) {
  if (x.norm() > d)
    x = d*(x/x.norm());
	
  return x;
}

float smoothstep(double x, double edge0, double edge1) {
  // Scale, bias and saturate x to 0..1 range
  x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
  // Evaluate polynomial
  return x * x * (3 - 2 * x);
}

double scaled_value(double x, double x_min, double x_max) {
    return (x - x_min) / (x_max - x_min);
}

VectorXd scaled_values(VectorXd const &x_vec) {

    double x_min = x_vec.minCoeff();
    double x_max = x_vec.maxCoeff();

    VectorXd y(x_vec.rows());

    for (int i = 0; i < x_vec.rows(); i++){
      y(i) = scaled_value(x_vec(i),x_min,x_max);
    }

    return y;

    //return x_vec.unaryExpr([this](double x) { return scaled_value(x,a,b); }).transpose();
}

#endif
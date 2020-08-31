#ifndef SMOOTHSTEP_H
#define SMOOTHSTEP_H

#include <algorithm>
#include <cmath>
#include "Eigen/Dense"
using namespace Eigen;

using namespace std;
typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;

float clamp(float x, float lowerlimit, float upperlimit) {
  if (x < lowerlimit)
    x = lowerlimit;
  if (x > upperlimit)
    x = upperlimit;
  return x;
}

Vxf clampmag(Vxf x, float d) {
  if (x.norm() > d)
    x = d*(x/x.norm());
	
  return x;
}

float smoothstep(float x, float edge0, float edge1) {
  // Scale, bias and saturate x to 0..1 range
  x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
  // Evaluate polynomial
  return x * x * (3 - 2 * x);
}

#endif
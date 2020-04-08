#ifndef _PTHRES_HPP_
#define _PTHRES_HPP_

#include "Eigen/Dense"
using namespace Eigen;

typedef Eigen::VectorXf Vxf;

inline Vxf pthres(Vxf V){
	V = ((V.array().abs() > 0.0).select(V.array(),0.0)).cast<float>();
	return V;
}

#endif

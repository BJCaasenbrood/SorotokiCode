#ifndef DELTA_H
#define DELTA_H

#include "Eigen/Dense"
using namespace Eigen;

typedef Eigen::VectorXf Vxf;

//---------------------------------------------------
//-------------------- delta operator for k-th entry
//---------------------------------------------------
Vxf delta(const int k, const int N)
{
	Vxf d(N);
	d.setZero();
	d(k) = 1.0;
	return d;
}

#endif
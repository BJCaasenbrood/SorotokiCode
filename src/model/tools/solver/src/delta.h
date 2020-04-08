#ifndef DELTA_H
#define DELTA_H

#include "Eigen/Dense"
using namespace Eigen;

typedef Eigen::VectorXf Vxf;

//---------------------------------------------------
//-------------------- delta operator for k-th entry
//---------------------------------------------------
Vxf isoSO3(const int k, const int N)
{
	Vxf d(N);
	d = Vxf::Zero(N);
	d(k) = 1.0;
	return d;
}

#endif
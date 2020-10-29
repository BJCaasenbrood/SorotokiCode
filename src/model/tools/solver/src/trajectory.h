#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#define _USE_MATH_DEFINES
#include <cmath>
#include "Eigen/Dense"
#include "smoothstep.h"
using namespace Eigen;
#define PI 3.1415926

typedef Eigen::Matrix<float, 6, 6> M6f;
typedef Eigen::Matrix<float, 4, 4> M4f;
typedef Eigen::Matrix<float, 3, 3> M3f;
typedef Eigen::Matrix<float, 7, 1> V7f;
typedef Eigen::Matrix<float, 6, 1> V6f;
typedef Eigen::Matrix<float, 4, 1> V4f;
typedef Eigen::Matrix<float, 3, 1> V3f;

typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;

#define OMEGA 0.3

float a = 0.01;
float b = 0.03;

//---------------------------------------------------
//----------------- isomorphism between R3 and so(3)
//---------------------------------------------------
void pathSolve(float t, V7f &yd)
{
	Vxf y0(7), yt(7);

	y0.setZero();
	yt.setZero();

	y0 = yd;

	float r;

	r = a*(1.0 + (1.0/8)*sin(2*OMEGA*t)*sin(2*OMEGA*t));

	//yt(4) = a*sin(OMEGA*t);
	//yt(5) = b*sin(OMEGA*t);
	//yt(6) = b*cos(OMEGA*t);

	yd.noalias() = y0 + yt;
}

#endif
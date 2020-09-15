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

#define OMEGA 0.75

float a = 0.3;
float b = 0.1;

//---------------------------------------------------
//----------------- isomorphism between R3 and so(3)
//---------------------------------------------------
void pathSolve(float t, V7f &yd)
{
	Vxf y0(7), yt(7);

	y0.setZero();
	yt.setZero();

	y0 = yd;

	yt(5) = a*sin(OMEGA*t);
	//yt(4) = b*sin(OMEGA*t + 1.5*PI);
	yt(6) = a*cos(OMEGA*t);

	yd.noalias() = y0 + yt;
}

#endif
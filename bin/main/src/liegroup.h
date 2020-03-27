#ifndef LIEGROUP_H
#define LIEGROUP_H

#include "Eigen/Dense"
using namespace Eigen;

typedef Eigen::Matrix<float, 6, 6> M6f;
typedef Eigen::Matrix<float, 3, 3> M3f;
typedef Eigen::Matrix<float, 7, 1> V7f;
typedef Eigen::Matrix<float, 6, 1> V6f;
typedef Eigen::Matrix<float, 4, 1> V4f;
typedef Eigen::Matrix<float, 3, 1> V3f;

typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;

//---------------------------------------------------
//----------------- isomorphism between R3 and so(3)
//---------------------------------------------------
M3f isoSO3(V3f x)
{
	M3f y;
	y << 0,-x(2),x(1),x(2),0,-x(0),-x(1),x(0),0;
	return y;
	
}

//---------------------------------------------------
//------------- quaternions to rotation matrix SO(3)
//---------------------------------------------------
M3f quat2rot(V4f q)
{
	float w,x,y,z;
	float Rxx,Rxy,Rxz;
	float Ryx,Ryy,Ryz;
	float Rzx,Rzy,Rzz;

	M3f R;

	w = q(0); 
	x = q(1); 
	y = q(2); 
	z = q(3);

	Rxx = 1 - 2*(y*y + z*z); 
	Rxy = 2*(x*y - z*w); 
	Rxz = 2*(x*z + y*w); 
	Ryx = 2*(x*y + z*w); 
	Ryy = 1 - 2*(x*x + z*z); 
	Ryz = 2*(y*z - x*w );
	Rzx = 2*(x*z - y*w ); 
	Rzy = 2*(y*z + x*w ); 
	Rzz = 1 - 2 *(x*x + y*y);

	R << Rxx, Rxy, Rxz, 
	     Ryx, Ryy, Ryz,
	     Rzx, Rzy, Rzz;

	return R;
}

//---------------------------------------------------
//-------------- adjoint action on Lie algebra se(3)
//---------------------------------------------------
M6f admap(V6f x)
{
	M6f ad;
	V3f W,U;
	W << x(0),x(1),x(2);
	U << x(3),x(4),x(5);

	M3f Wh = isoSO3(W);
	M3f Uh = isoSO3(U);

	ad.block(0,0,3,3) = Wh;
	ad.block(3,3,3,3) = Wh;
	ad.block(3,0,3,3) = Uh;
	ad.block(0,3,3,3) = M3f::Zero(3,3);

	return ad;

}

//---------------------------------------------------
//------------------ adjoint map on lie group SE(3)
//---------------------------------------------------
M6f Admap(V7f x)
{
	M6f Ad;
	M3f R;
	V4f quat;
	V3f r;
	
	quat << x(0),x(1),x(2),x(3);
	r << x(4),x(5),x(6);

	R = quat2rot(quat);
	Ad.block(0,0,3,3) = R;
	Ad.block(3,3,3,3) = R;
	Ad.block(3,0,3,3) = isoSO3(r)*R;
	Ad.block(0,3,3,3) = M3f::Zero(3,3);

	return Ad;
}

#endif
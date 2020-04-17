#ifndef LIEGROUP_H
#define LIEGROUP_H

#include "Eigen/Dense"
using namespace Eigen;

typedef Eigen::Matrix<float, 6, 6> M6f;
typedef Eigen::Matrix<float, 4, 4> M4f;
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
void isoSO3(V3f x, M3f &y)
{
	y << 0,-x(2),x(1),x(2),0,-x(0),-x(1),x(0),0;
}

//---------------------------------------------------
//------------- quaternions to rotation matrix SO(3)
//---------------------------------------------------
void quat2rot(V4f q, M3f &R)
{
	float w,x,y,z;
	float Rxx,Rxy,Rxz;
	float Ryx,Ryy,Ryz;
	float Rzx,Rzy,Rzz;

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

}

//---------------------------------------------------
//-------------- adjoint action on Lie algebra se(3)
//---------------------------------------------------
void admap(V6f x, M6f &ad)
{
	V3f W;
	M3f Wh;

	W << x(0),x(1),x(2);
	isoSO3(W,Wh);

	ad.block(0,0,3,3).noalias() = Wh;
	ad.block(3,3,3,3).noalias() = Wh;

	W << x(3),x(4),x(5);
	isoSO3(W,Wh);
	ad.block(3,0,3,3).noalias() = Wh;
	ad.block(0,3,3,3).noalias() = M3f::Zero(3,3);
}

//---------------------------------------------------
//------------------ adjoint map on lie group SE(3)
//---------------------------------------------------
M6f Admap(V7f x)
{
	M6f Ad;
	M3f R,S;
	V4f w;
	V3f v;
	
	w << x(0),x(1),x(2),x(3);
	v << x(4),x(5),x(6);

	quat2rot(w,R);
	isoSO3(v,S);

	Ad.block(0,0,3,3).noalias() = R;
	Ad.block(3,3,3,3).noalias() = R;
	Ad.block(0,3,3,3).noalias() = M3f::Zero(3,3);
	(Ad.block(3,0,3,3)).noalias() = S*R;

	return Ad;
}

//---------------------------------------------------
//------------------ adjoint map on lie group SE(3)
//---------------------------------------------------
M6f AdmapInv(V7f x)
{
	M6f Ad;
	M3f R, S;
	V4f w;
	V3f v;
	
	w << x(0),x(1),x(2),x(3);
	v << x(4),x(5),x(6);

	quat2rot(w,R);
	R.transposeInPlace();

	isoSO3(v,S);

	Ad.block(0,0,3,3).noalias() = R;
	Ad.block(3,3,3,3).noalias() = R;
	Ad.block(0,3,3,3).noalias() = M3f::Zero(3,3);
	(Ad.block(3,0,3,3)).noalias() = R*S.transpose();

	return Ad;
}

//---------------------------------------------------
//--------------- curvature-twist to matrix operator
//---------------------------------------------------
void strainMapping(V3f k, M4f &K){

	float k1 = k(0);
	float k2 = k(1);
	float k3 = k(2);

	K <<  0, -k1, -k2, -k3, 
	     k1,   0, -k3,  k2, 
         k2,  k3,   0, -k1, 
         k3, -k2,  k1,   0;

}

#endif
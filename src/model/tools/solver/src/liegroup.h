#ifndef LIEGROUP_H
#define LIEGROUP_H

#define _USE_MATH_DEFINES
#include <cmath>
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
void SE3toR6(V7f g, V6f &t)
{

	float w,x,y,z;
	float roll,pitch,yaw;
	w = g(0); 
	x = g(1); 
	y = g(2); 
	z = g(3);

    float sinr_cosp = 2*(w*x + y*z);
    float cosr_cosp = 1- 2*(x*x + y*y);
    roll = std::atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    float sinp = 2 * (w*y - z*x);
    if (std::abs(sinp) >= 1){
        pitch = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    }
    else{
        pitch= std::asin(sinp);
    }

    // yaw (z-axis rotation)
    float siny_cosp = 2 * (w*z + x*y);
    float cosy_cosp = 1 - 2*(y*y + z*z);
    yaw = std::atan2(siny_cosp, cosy_cosp);

    t << roll,pitch,yaw,g(4),g(5),g(6);
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
void Admap(V7f x, M6f &Ad)
{
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
}

//---------------------------------------------------
//------------------ adjoint map on lie group SE(3)
//---------------------------------------------------
void AdmapInv(V7f x, M6f &Ad)
{
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
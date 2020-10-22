#ifndef LIEGROUP_H
#define LIEGROUP_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
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
//-------------------------------------------- SE(3)
//---------------------------------------------------
void SE3(V7f g, M4f &X)
{

	M3f R;
	V4f q;

	q << g(0),g(1),g(2),g(3);
	quat2rot(q,R);

	X(3,0) = 0.0;
	X(3,1) = 0.0;
	X(3,2) = 0.0;
	X(3,3) = 1.0;
	X.block(0,0,3,3).noalias() = R;
	X.block(0,3,3,1).noalias() = g.block(4,0,3,1);
}

//---------------------------------------------------
//-------------------------------------------- SE(3)
//---------------------------------------------------
void SE3Inv(V7f g, M4f &X)
{

	M3f R;
	V4f q;

	q << g(0),g(1),g(2),g(3);
	quat2rot(q,R);

	X(3,0) = 0.0;
	X(3,1) = 0.0;
	X(3,2) = 0.0;
	X(3,3) = 1.0;
	X.block(0,0,3,3).noalias() = R.transpose();
	X.block(0,3,3,1).noalias() = -R.transpose()*g.block(4,0,3,1);
}

//---------------------------------------------------
//----------------- isomorphism between R3 and so(3)
//---------------------------------------------------
void isoSO3(V3f x, M3f &y)
{
	y << 0.0,-x(2),x(1),x(2),0,-x(0),-x(1),x(0),0.0;
}

//---------------------------------------------------
//----------------- isomorphism between so(3) and R3
//---------------------------------------------------
void isoSO3_(M3f x, V3f &y)
{
	y << x(2,1),x(0,2),x(1,0);
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
        pitch = std::copysign(M_PI / 2, sinp);
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
//-------------- adjoint action on Lie algebra se(3)
//---------------------------------------------------
void coadmap(V6f x, M6f &coad)
{
	V3f W;
	M3f Wh;

	W << x(0),x(1),x(2);
	isoSO3(W,Wh);

	coad.block(0,0,3,3).noalias() = Wh;
	coad.block(3,3,3,3).noalias() = Wh;

	W << x(3),x(4),x(5);
	isoSO3(W,Wh);
	coad.block(0,3,3,3).noalias() = Wh;
	coad.block(3,0,3,3).noalias() = M3f::Zero(3,3);
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

//---------------------------------------------------
//------------------------- exponential map of SO(3)
//---------------------------------------------------
void expmapSO3(V3f x, M3f &Y){

	M3f I, X;

	I.noalias() = M3f::Identity(3,3);
	isoSO3(x,X);

	float t = x.norm();
	float a = 1.0;
	float b = 1.0;

	if(abs(t) >= 1e-6){
		a = sin(t)/t;
		b = (1.0-cos(t))/(t*t);
	}

	Y.noalias() = I + a*X + b*X*X;

}

//---------------------------------------------------
//------------------------- logarithmic map of SO(3)
//---------------------------------------------------
void logmapSO3(M3f X, M3f &Y){


	float theta = acos(0.5*(X.trace()-1.0));
	float a = 0.5*theta/sin(theta);

	if(!isnan(a)){
		Y.noalias() = a*(X - X.transpose());
	}
	else{
		Y.setZero();
	}

}

//---------------------------------------------------
//-------------------------- tangential map of SO(3)
//---------------------------------------------------
void tmapSO3(V3f x, M3f &Y){

	M3f I, X;

	I.noalias() = M3f::Identity(3,3);
	isoSO3(x,X);

	float t = x.norm();
	float a = 1.0;
	float b = 1.0;

	if(abs(t) >= 1e-6){
		a = sin(t)/t;
		b = (1.0-cos(t))/(t*t);
	}

	if(abs(t) >= 1e-6){
		Y.noalias() = I - b*X + (1/(t*t))*(1-a)*X*X;
	}
	else{
		Y.noalias() = I;
	}
	

}

//---------------------------------------------------
//------------------ inverse tangential map of SO(3)
//---------------------------------------------------
void tmapInvSO3(V3f x, M3f &Y){

	M3f I, X;

	I.noalias() = M3f::Identity(3,3);
	isoSO3(x,X);

	float t = x.norm();
	float a = 1.0;
	float b = 1.0;

	if(abs(t) >= 1e-6){
		a = sin(t)/t;
		b = 2.0*(1.0-cos(t))/(t*t);
	}

	if(abs(t) >= 1e-6){
		Y.noalias() = I + 0.5*X + (1/(t*t))*(1-a/b)*X*X;
	}
	else{
		Y.noalias() = I;
	}

}


//---------------------------------------------------
//------------------------- exponential map of SE(3)
//---------------------------------------------------
void expmapSE3(V6f x, M4f &Y){

	V3f W;
	M3f Wh;

	Y(3,0) = 0.0;
	Y(3,1) = 0.0;
	Y(3,2) = 0.0;
	Y(3,3) = 1.0;

	W << x(0),x(1),x(2);
	expmapSO3(W,Wh);

	Y.block(0,0,3,3).noalias() = Wh;

	tmapSO3(W,Wh);
	W << x(3),x(4),x(5);

	Y.block(0,3,3,1).noalias() = Wh.transpose()*W;

}

//---------------------------------------------------
//------------------------- logarithmic map of SE(3)
//---------------------------------------------------
void logmapSE3(M4f X, V6f &Y){

	V3f W;
	M3f R,Wh;

	R << X(0,0),X(0,1),X(0,2),
		 X(1,0),X(1,1),X(1,2),
		 X(2,0),X(2,1),X(2,2);

	logmapSO3(R,Wh);
	isoSO3_(Wh,W);

	// IS THIS CORRECT?
	Y.block(0,0,3,1).noalias() = W;

	tmapInvSO3(W,Wh);

	W << X(0,3),X(1,3),X(2,3);

	// IS THIS CORRECT?
	Y.block(3,0,3,1).noalias() = Wh.transpose()*W;

}

//---------------------------------------------------
//------------------ omega-U-tangential map of SO(3)
//---------------------------------------------------
void tou_plus_SE3(V3f xu, V3f xo, M3f &T){

	M3f A,B,C;

	isoSO3(xu,A);
	isoSO3(xo,B);

	C.noalias() = A*B + B*A;

	float t = xo.norm();
	float a = 1.0;
	float b = 1.0;

	if(abs(t) >= 1e-6){
		a = sin(t)/t;
		b = 2.0*(1.0-cos(t))/(t*t);
	}

	if(abs(t) <= 1e-6){
		T.noalias() = -0.5*A;
	}
	else{
		T.noalias() = -0.5*b*A + (1/(t*t))*(1-a)*C 
			+ (1/(t*t))*B.transpose()*A*((b-a)*B 
			+ (0.5*b - (1/(t*t))*(3*(1-a)))*B*B);
	}

}

//---------------------------------------------------
//-------------------------- tangential map of SO(3)
//---------------------------------------------------
void tmapSE3(V6f x, M6f &T){

	M3f To,Tou;
	V3f V,W;

	W << x(0),x(1),x(2);
	V << x(3),x(4),x(5);
	tmapSO3(W,To);
	tou_plus_SE3(V,W,Tou);

	T.block(0,0,3,3).noalias() = To;
	T.block(3,3,3,3).noalias() = To;
	T.block(3,0,3,3).noalias() = M3f::Zero(3,3);
	T.block(0,3,3,3).noalias() = Tou;

}


#endif
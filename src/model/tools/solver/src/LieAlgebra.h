#ifndef LIEGROUP_H
#define LIEGROUP_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
using namespace Eigen;

typedef Eigen::Matrix<double, 6, 6> M6d;
typedef Eigen::Matrix<double, 4, 4> M4d;
typedef Eigen::Matrix<double, 3, 3> M3d;
typedef Eigen::Vector<double, 13> V13d;
typedef Eigen::Vector<double, 7> V7d;
typedef Eigen::Vector<double, 6> V6d;
typedef Eigen::Vector<double, 4> V4d;
typedef Eigen::Vector<double, 3> V3d;


//---------------------------------------------------
//------------- quaternions to rotation matrix SO(3)
//---------------------------------------------------
void quat2rot(V4d q, M3d &R)
{
	double w,x,y,z;
	double Rxx,Rxy,Rxz;
	double Ryx,Ryy,Ryz;
	double Rzx,Rzy,Rzz;

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
void SE3(V7d g, M4d &X)
{

	M3d R;
	V4d q;

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
void SE3Inv(V7d g, M4d &X)
{

	M3d R;
	V4d q;

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
void isoSO3(V3d x, M3d &y)
{
	y << 0.0,-x(2),x(1),x(2),0,-x(0),-x(1),x(0),0.0;
}

//---------------------------------------------------
//----------------- isomorphism between so(3) and R3
//---------------------------------------------------
void isoSO3_(M3d x, V3d &y)
{
	y << x(2,1),x(0,2),x(1,0);
}

//---------------------------------------------------
//---- conversion from SE(3) to (xyz)+(Euler-angles)
//---------------------------------------------------
void SE3toR6(V7d g, V6d &t)
{

	double w,x,y,z;
	double roll,pitch,yaw;
	w = g(0); 
	x = g(1); 
	y = g(2); 
	z = g(3);

    double sinr_cosp = 2*(w*x + y*z);
    double cosr_cosp = 1- 2*(x*x + y*y);
    roll = std::atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    double sinp = 2 * (w*y - z*x);
    if (std::abs(sinp) >= 1){
        pitch = std::copysign(M_PI / 2, sinp);
    }
    else{
        pitch= std::asin(sinp);
    }

    // yaw (z-axis rotation)
    double siny_cosp = 2 * (w*z + x*y);
    double cosy_cosp = 1 - 2*(y*y + z*z);
    yaw = std::atan2(siny_cosp, cosy_cosp);

    t << roll,pitch,yaw,g(4),g(5),g(6);
}

//---------------------------------------------------
//---- conversion from SE(3) to (xyz)
//---------------------------------------------------
void SE3pos(V7d g, V6d &t)
{
    t << 0,0,0,g(4),g(5),g(6);
}


//---------------------------------------------------
//-------------- adjoint action on Lie algebra se(3)
//---------------------------------------------------
void admap(V6d x, M6d &ad)
{
	V3d W;
	M3d Wh;

	W << x(0),x(1),x(2);
	isoSO3(W,Wh);

	ad.block(0,0,3,3).noalias() = Wh;
	ad.block(3,3,3,3).noalias() = Wh;

	W << x(3),x(4),x(5);
	isoSO3(W,Wh);
	ad.block(3,0,3,3).noalias() = Wh;
	ad.block(0,3,3,3).noalias() = M3d::Zero(3,3);
}

//---------------------------------------------------
//-------------- adjoint action on Lie algebra se(3)
//---------------------------------------------------
void coadmap(V6d x, M6d &coad)
{
	V3d W;
	M3d Wh;

	W << x(0),x(1),x(2);
	isoSO3(W,Wh);

	coad.block(0,0,3,3).noalias() = Wh;
	coad.block(3,3,3,3).noalias() = Wh;

	W << x(3),x(4),x(5);
	isoSO3(W,Wh);
	coad.block(0,3,3,3).noalias() = Wh;
	coad.block(3,0,3,3).noalias() = M3d::Zero(3,3);
}

//---------------------------------------------------
//------------------ adjoint map on lie group SE(3)
//---------------------------------------------------
void Admap(V7d x, M6d &Ad)
{
	M3d R,S;
	V4d w;
	V3d v;
	
	w << x(0),x(1),x(2),x(3);
	v << x(4),x(5),x(6);

	quat2rot(w,R);
	isoSO3(v,S);

	Ad.block(0,0,3,3).noalias() = R;
	Ad.block(3,3,3,3).noalias() = R;
	Ad.block(0,3,3,3).noalias() = M3d::Zero(3,3);
	(Ad.block(3,0,3,3)).noalias() = S*R;
}

//---------------------------------------------------
//------------------ adjoint map on lie group SE(3)
//---------------------------------------------------
void AdmapInv(V7d x, M6d &Ad)
{
	M3d R, S;
	V4d w;
	V3d v;
	
	w << x(0),x(1),x(2),x(3);
	v << x(4),x(5),x(6);

	quat2rot(w,R);
	R.transposeInPlace();

	isoSO3(v,S);

	Ad.block(0,0,3,3).noalias() = R;
	Ad.block(3,3,3,3).noalias() = R;
	Ad.block(0,3,3,3).noalias() = M3d::Zero(3,3);
	(Ad.block(3,0,3,3)).noalias() = R*S.transpose();
}

//---------------------------------------------------
//--------------- curvature-twist to matrix operator
//---------------------------------------------------
void strainMapping(V3d k, M4d &K){

	double k1 = k(0);
	double k2 = k(1);
	double k3 = k(2);

	K <<  0, -k1, -k2, -k3, 
	     k1,   0, -k3,  k2, 
         k2,  k3,   0, -k1, 
         k3, -k2,  k1,   0;

}

//---------------------------------------------------
//------------------------- exponential map of SO(3)
//---------------------------------------------------
void expmapSO3(V3d x, M3d &Y){

	M3d X;

	//I.noalias() = M3f::Identity(3,3);
	isoSO3(x,X);

	double t = x.norm();
	double a = 1.0;
	double b = 1.0;

	if(abs(t) >= 1e-6){
		a = sin(t)/t;
		b = (1.0-cos(t))/(t*t);
	}

	Y.noalias() = M3d::Identity(3,3) + a*X + b*X*X;

}

//---------------------------------------------------
//------------------------- logarithmic map of SO(3)
//---------------------------------------------------
void logmapSO3(M3d X, M3d &Y){


	double theta = acos(0.5*(X.trace()-1.0));
	double a = 0.5*theta/sin(theta);

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
void tmapSO3(V3d x, M3d &Y){

	M3d X;

	//I.noalias() = M3f::Identity(3,3);
	isoSO3(x,X);

	double t = x.norm();
	double a = 1.0;
	double b = 1.0;

	if(abs(t) >= 1e-6){
		a = sin(t)/t;
		b = (1.0-cos(t))/(t*t);
	}

	if(abs(t) >= 1e-6){
		Y.noalias() = M3d::Identity(3,3) - b*X + (1/(t*t))*(1-a)*X*X;
	}
	else{
		Y.noalias() = M3d::Identity(3,3);
	}
	

}

//---------------------------------------------------
//------------------ inverse tangential map of SO(3)
//---------------------------------------------------
void tmapInvSO3(V3d x, M3d &Y){

	M3d X;

	//I.noalias() = M3f::Identity(3,3);
	isoSO3(x,X);

	double t = x.norm();
	double a = 1.0;
	double b = 1.0;

	if(abs(t) >= 1e-6){
		a = sin(t)/t;
		b = 2.0*(1.0-cos(t))/(t*t);
	}

	if(abs(t) >= 1e-6){
		Y.noalias() = M3d::Identity(3,3) + 0.5*X + (1/(t*t))*(1-a/b)*X*X;
	}
	else{
		Y.noalias() = M3d::Identity(3,3);
	}

}


//---------------------------------------------------
//------------------------- exponential map of SE(3)
//---------------------------------------------------
void expmapSE3(V6d x, M4d &Y){

	V3d W;
	M3d Wh;

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
void logmapSE3(M4d X, V6d &Y){

	V3d W;
	M3d R,Wh;

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
void tou_plus_SE3(V3d xu, V3d xo, M3d &T){

	M3d A,B,C;

	isoSO3(xu,A);
	isoSO3(xo,B);

	C.noalias() = A*B + B*A;

	double t = xo.norm();
	double a = 1.0;
	double b = 1.0;

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
void tmapSE3(V6d x, M6d &T){

	M3d To,Tou;
	V3d V,W;

	W << x(0),x(1),x(2);
	V << x(3),x(4),x(5);
	tmapSO3(W,To);
	tou_plus_SE3(V,W,Tou);

	T.block(0,0,3,3).noalias() = To;
	T.block(3,3,3,3).noalias() = To;
	T.block(3,0,3,3).noalias() = M3d::Zero(3,3);
	T.block(0,3,3,3).noalias() = Tou;

}


#endif
#ifndef MODEL_H
#define MODEL_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include "liegroup.h"
#include "shapesx.cpp"
#include "qprog.hpp"
#include "pinv.cpp"
#include "tictoc.h"
#include "smoothstep.h"
#include "ConfigFile.cpp"
#include "Chameleon.cpp"

//#define FULL_CONTROLLER
//#define CONSTRAINED_CONTROLLER
//#define ENERGY_PROJECTION_CONTROLLER
#define SOLVER_OUTPUT
//#define WRITE_OUTPUT
#define TICTOC
//#define QUASINETWON
#define JACOBIAN

#define PRECISION 5
#define PI 3.1415926

typedef Eigen::Array<int, Dynamic, 1> Vxi;
typedef Eigen::Array<int, 6, 1> V6i;
typedef Eigen::Matrix<float, 6, 6> M6f;
typedef Eigen::Matrix<float, 4, 4> M4f;
typedef Eigen::Matrix<float, 3, 3> M3f;
typedef Eigen::Matrix<float, 7, 1> V7f;
typedef Eigen::Matrix<float, 13, 1> V13f;
typedef Eigen::Matrix<float, 6, 1> V6f;
typedef Eigen::Matrix<float, 4, 1> V4f;
typedef Eigen::Matrix<float, 3, 1> V3f;
typedef Eigen::Matrix<float, 2, 1> V2f;
typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;
typedef Eigen::VectorXd Vxd;
typedef Eigen::MatrixXd Mxd;

typedef Eigen::Matrix<float, 19, 1> Vff;

class Model
{
  public:

  	bool ENERGY_CONTROLLER;
  	bool WRITE_OUTPUT;

  	int E1,E2,E3;
  	int K1,K2,K3; 

  	float E11, E22, E33;
  	float G11, G22, G33;

  	int NMODE; 
  	int NDISC;
	float SDOMAIN; 
	float TDOMAIN; 
	int SPACESTEP;
	float TIMESTEP; 
	int INTSTEP;

	float ATOL;
	float RTOL; 
	int MAX_IMPL;
	int MAX_ITER; 
	float SPEEDUP; 

	float PRS_AREA; 
	float GRAVITY; 
	float RADIUS; 
	float RHO; 
	float EMOD; 
	float NU; 
	float MU; 

  	float KP, KD, KE;
  	float Ec;
  	int NDof, NState;
  	Shapes Phi;

  	Mxf Hess;
  	V2f Hamil; // H:=[T,V]

  	Mxf Ba,Bc;
  	Mxf Sa,Sc;
  	M6f Mtt, Ktt, Dtt;
  	Mxf Mee, Kee, Dee;
  	Vxf q, dq, ddq;
  	Vxf Qa, Qv, Qu, Qd;
  	Vxf qd;
  	V6f Xi0;

  	V7f g;
  	V6f xi,dxi,ddxi,eta,deta,lam;
	V4f quat;
	V3f xiK, xiE;
	M3f R; 
	M4f A;

	V6f gvec;
	Vxf z0, tau, u;

	Model(const char* str);

	void output(ofstream &file, float t, Vxf x);
	void read(const char* str, Vxf &x);
	void cleanup();

	//void operationalSpaceDynamics(Mxf &J, Mxf &Jt, Vxf &dq, 
	//Mxf &M, Vxf &C, Vxf &G, Mxf &Mx, Vxf &Cx, Vxf &Gx);
	//void dynamicProjector(Mxf &J, Mxf &M, Mxf &S, Mxf &P);
	//void controllerWrench(float t, Mxf &J, Vxf &f);
	void controllerPassive(float t, Mxf &M, Vxf &Qv,
 	Vxf &Qa, Mxf &J, Vxf &f);

	/*
	Vxf solve();
	Vxf implicit_solve();
	Vxf simulate();
	*/
	
	Vxf implicit_simulate();

	void inverseDynamics(Vxf v, Vxf dv, Vxf ddv, Vxf &Q);
	void buildJacobian(float se, Mxf &J, Mxf &Jt);
	void buildInertia(Vxf x, Mxf &M);

	void kinematicODE(float t, Vxf x, Vxf &dx);
	void dynamicODE(float t, Vxf x, Vxf &dx);

	void forwardODE(float s, Vff x, Vff &dx);
	void backwardODE(float s, Vxf x, Vxf &dx);
	void jacobiODE(float s, V13f x, V13f &dx, Mxf &dJ, Mxf &dJt);
	void inertiaODE(float s, V7f x, Mxf J, 
	V7f &dx, Mxf &dJ, Mxf &dM);
	void systemMatODE(float s, 
	Mxf &K, Mxf &M, Mxf &D);

	Mxf tableConstraints(Vxi table, bool set = true);
	void hessianInverse(float dt, Vxf R, Vxf &dx);

	void buildInertiaTensor();
	void buildStiffnessTensor();
	void buildDampingTensor();
	void buildGlobalSystem();

	Mxf pressureMapping();
};

#endif
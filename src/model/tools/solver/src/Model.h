#ifndef MODEL_H
#define MODEL_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include "liegroup.h"
#include "shapes.cpp"
#include "qprog.hpp"
#include "tictoc.h"
#include "smoothstep.h"

//#define FULL_CONTROLLER
//#define CONSTRAINED_CONTROLLER
#define WRITE_OUTPUT
#define TICTOC
//#define QUASINEWTON
//#define JACOBIAN

#define PRECISION 5

#define NMODE 3
#define SDOMAIN 1
#define TDOMAIN 10
#define SPACESTEP 10
#define TIMESTEP 0.0333
#define INTSTEP 900

#define ATOL 1e-4
#define RTOL 1e-3
#define MAX_ITER 1e5
#define SPEEDUP 1.0

#define PRS_AREA 1e-5
#define GRAVITY 9.81
#define PI 3.1415926
#define RADIUS 0.01
#define RHO 0.01
#define EMOD 4e2
#define NU 0.4
#define MU 0.2
//#define MU 0.01 // (static solver)

typedef Eigen::Array<int, 6, 1> V6i;
typedef Eigen::Matrix<float, 6, 6> M6f;
typedef Eigen::Matrix<float, 4, 4> M4f;
typedef Eigen::Matrix<float, 3, 3> M3f;
typedef Eigen::Matrix<float, 7, 1> V7f;
typedef Eigen::Matrix<float, 13, 1> V13f;
typedef Eigen::Matrix<float, 6, 1> V6f;
typedef Eigen::Matrix<float, 4, 1> V4f;
typedef Eigen::Matrix<float, 3, 1> V3f;
typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;
typedef Eigen::VectorXd Vxd;
typedef Eigen::MatrixXd Mxd;

typedef Eigen::Matrix<float, 19, 1> Vff;

class Model
{
  public:

  	int NDof, NState;
  	Shapes Phi;

  	Mxf Hess;

  	Mxf Ba;
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
	Vxf tau;

	Model(V6i table);

	void output(const char* str, float t, Vxf x);
	void read(const char* str, Vxf &x);
	void cleanup();

	void realtimeController(float t, Mxf &A1, Mxf &A2, Vxf Qdes, Vxf &X);

	Vxf solve();
	Vxf implicit_solve();
	Vxf simulate();
	Vxf implicit_simulate();

	void inverseDynamics(Vxf v, Vxf dv, Vxf ddv, Vxf &Q);
	void buildJacobian(float se, Mxf &J);
	void buildInertia(Vxf x, Mxf &M);

	void kinematicODE(float t, Vxf x, Vxf &dx);
	void dynamicODE(float t, Vxf x, Vxf &dx);

	void forwardODE(float s, Vff x, Vff &dx);
	void backwardODE(float s, Vxf x, Vxf &dx);
	void jacobiODE(float s, V13f x, V13f &dx, Mxf &dJ);
	void inertiaODE(float s, V7f x, Mxf J, 
	V7f &dx, Mxf &dJ, Mxf &dM);
	void systemMatODE(float s, 
	Mxf &K, Mxf &M, Mxf &D);

	Mxf tableConstraints(V6i table);
	void hessianInverse(float dt, Vxf R, Vxf &dx);

	void buildInertiaTensor();
	void buildStiffnessTensor();
	void buildDampingTensor();
	void buildGlobalSystem();

	Mxf pressureMapping();
};

#endif
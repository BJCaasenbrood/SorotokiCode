#ifndef MODEL_H
#define MODEL_H

#include <cstdio>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include <cmath>
#include <vector> 
#include "liegroup.h"
#include "shapes.cpp"
#include "qprog.hpp"
#include "pthres.hpp"
#include "tictoc.h"
#include "smoothstep.h"

#define WRITE_OUTPUT
#define TICTOC
//#define QUASINEWTON

#define PRECISION 5

#define NMODE 3
#define SDOMAIN 1
#define TDOMAIN 10
#define SPACESTEP 15
#define TIMESTEP 0.05
#define INTSTEP 900

#define ATOL 1e-6
#define RTOL 1e-4
#define MAX_ITER 1e5
#define SPEEDUP 1.0

#define PRS_AREA 1e-5
#define GRAVITY 0.00
#define PI 3.1415926
#define RADIUS 0.01
#define RHO 0.005
#define EMOD 400
#define NU 0.4
#define MU 0.3
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
	void read();
	void cleanup();

	void realtimeController(float t, Mxf J1, Mxf J2);

	Vxf solve();
	Vxf implicit_solve();
	Vxf simulate();
	Vxf implicit_simulate();
	void fmincon(Mxf &A, Vxf &b, Vxf &x);

	Vxf inverseDynamics(Vxf v, Vxf dv, Vxf ddv);
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
	Mxf hessianInverse(float dt);

	void buildInertiaTensor();
	void buildStiffnessTensor();
	void buildDampingTensor();
	void buildGlobalSystem();
	
	M4f strainMapping(V3f k);
	Mxf pressureMapping();
};

#endif
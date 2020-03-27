#ifndef MODEL_H
#define MODEL_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <cmath>
#include "liegroup.h"
#include "shapes.cpp"

#define SDOMAIN 1
#define TDOMAIN 10
#define SPACESTEP 100
#define TIMESTEP 5
#define INTSTEP 9000

#define PRECISION 5
#define TOLERANCE 1e-5
#define MAX_ITER 1e5
#define ADAPTIVE false
#define ADAP_KP 0.25
#define ADAP_KI 0.175
#define ADAP_MAX 0.001
#define ADAP_MIN 1.0

#define PRS_AREA 1e-5
#define GRAVITY 9.81
#define PI 3.1415926
#define RADIUS 0.01
#define RHO 0.01
#define EMOD 400
#define NU 0.4
#define MU 0.5

typedef Eigen::Array<int, 6, 1> V6i;
typedef Eigen::Matrix<float, 6, 6> M6f;
typedef Eigen::Matrix<float, 4, 4> M4f;
typedef Eigen::Matrix<float, 3, 3> M3f;
typedef Eigen::Matrix<float, 7, 1> V7f;
typedef Eigen::Matrix<float, 6, 1> V6f;
typedef Eigen::Matrix<float, 4, 1> V4f;
typedef Eigen::Matrix<float, 3, 1> V3f;
typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;

typedef Eigen::Matrix<float, 19, 1> Vff;

class Model
{
  public:

  	int NMode, NDof, NState;
  	Shapes Phi;

  	Mxf Hess;

  	Mxf Ba;
  	M6f Mtt, Ktt, Dtt;
  	Mxf Mee, Kee, Dee;
  	Vxf q, dq, ddq;
  	V6f Xi0;

	V3f gvec;
	V3f P1,P2;

	Model(V6i table, int nmode);

	Vxf solve();
	Vxf implicit_solve();
	Vxf inverseDynamics(Vxf v, Vxf dv, Vxf ddv);
	Mxf buildJacobian(float se);

	void kinematicODE(float t, Vxf x, Vxf &dx);
	void forwardODE(float s, Vff x, Vff &dx);
	void backwardODE(float s, Vxf x, Vxf &dx);
	void jacobiODE(float s, V7f x, V7f &dx, Mxf &dJ);

	Mxf tableConstraints(V6i table);
	void buildStiffnessTensor();
	void buildInertiaTensor();
	void buildDampingTensor();
	void buildGlobalSystem();

	void buildJacobian(float s, Vff x);
	
	M4f strainMapping(V3f k);
	Mxf pressureMapping();
};

#endif
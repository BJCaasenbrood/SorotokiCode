#ifndef MODEL_PH_H
#define MODEL_PH_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <unistd.h>

#include "lieAlgebra.h"
#include "table.h"
#include "tictoc.h"
#include "smoothstep.h"
#include "trajectory.h"
#include "readFile.h"
#include "Config/ConfigFile.cpp"
#include "Config/Chameleon.cpp"

#ifdef DISCONTINIOUS
#include "shapesx.cpp"
#else
#include "shapes.cpp"
#endif

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
typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;
typedef Eigen::VectorXd Vxd;
typedef Eigen::MatrixXd Mxd;

typedef Eigen::Matrix<float, 19, 1> Vff;

typedef Eigen::CompleteOrthogonalDecomposition<Matrix<float, Dynamic, Dynamic>> Cxd;

class Model
{
  public:

  	bool ENERGY_CONTROLLER;
  	bool KINEMATIC_CONTROLLER;
  	bool WRITE_OUTPUT;
  	bool POINT_INPUT;

  	float A11, J11, J22, J33;
 	float E11, E22, E33;
  	float G11, G22, G33;

  	int NMODE; 
  	int NDISC;
	float SDOMAIN; 
	float TDOMAIN; 
	int SPACESTEP;
	float TIMESTEP; 
	int INTSTEP;

	float ATOL, RTOL, KTOL; 
	int MAX_IMPL, MAX_ITER, MAX_IK;
	float LAMBDA;
	float SPEEDUP; 

	float PRS_AREA; 
	float GRAVITY; 
	float RADIUS; 
	float RHO; 
	float EMOD; 
	float NU; 
	float MU; 

	float ADAMP;

  	float KP, KD;
  	float Hstar;
  	int NDof, NState;
  	Shapes Phi;

  	Mxf Hess;

  	Mxf Ba, Bc;
  	Mxf Sa, Sc;
  	M6f Mtt, Ktt, Dtt;
  	Mxf Mtee, Mee, Cee, Kee, Dee;
  	Mxf Mbee, Cbee;
  	Vxf Gee;
  	Vxf q, dq, ddq;
  	Vxf dq_;
  	Vxf Qa, Qv, Qu, Qd;
  	Vxf qd;
  	V6f Xi0;

  	V7f g, gd;
  	V6f xi,dxi,ddxi,eta,deta,lam;
	V4f quat;
	V3f xiK, xiE;
	M3f R; 
	M4f A;

	Vxf gvec;
	Vxf z0, tau, u;

	Model(const char* str);

	void output(ofstream &file, float t, Vxf x);
	void read(const char* str, Vxf &x);
	void cleanup();

	void controllerPassive(float t, Vxf &Hq, Vxf &Hp, Vxf &f);

	void inverse_kinematics(float t = 1.0);
	Vxf simulate();
	Vxf implicit_simulate();

	void buildJacobian(float se, Mxf &J, Mxf &Jt);
	void buildLagrange(Vxf v, Vxf dv, Mxf &M, Mxf &C,
	Vxf &dG, Mxf &Mt, Mxf &Me, Mxf &Ce);
	
	void dynamicODE(float t, Vxf x, Vxf &dx);
	void jacobiODE(float s, V13f x, V13f &dx, Mxf &dJ, Mxf &dJt);
	void systemMatODE(float s,Mxf &K, Mxf &M, Mxf &D);
	void lagrangianODE(float s, V13f x, Mxf J, Mxf Jt,
	V13f &dx, Mxf &dJ, Mxf &dJt, Mxf &dM, Mxf &dC,
	Vxf &dG, Mxf &dMt, Mxf &dMe, Mxf &dCe);

	void hessianInverse(float alpha, Vxf R, Vxf &dx);

	void buildInertiaTensor();
	void buildStiffnessTensor();
	void buildDampingTensor();
	void buildGlobalSystem();
	void buildNonlinearElastic(Vxf x, Vxf &N);

	Mxf pressureMapping();

	float Hamiltonian(Vxf x1, Vxf x2);
	void PoissonBracket(Vxf &Hq, Vxf &Hp, Mxf &A);
};

#endif
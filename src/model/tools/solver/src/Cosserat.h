#ifndef COSSERAT_H
#define COSSERAT_H

#include <Eigen/Dense>
#include "LieAlgebra.h"
#include "PropertyTensor.h"
#include "ConfigFile.cpp"
#include "Chameleon.cpp"
#include "table.h"
#include "readFile.h"
#include "DEBUGGER.h"

#define DISCONTINIOUS

#ifdef DISCONTINIOUS
#include "ShapesX.cpp"
#else
#include "Shapes.cpp"
#endif

using namespace std;
using namespace Eigen;

//template <short Na, short Nm >
class Cosserat 
{
	public:

	int Na;		// number of active dof
	int Nm;		// number of modes
	int Nd;		// number of discontinuities

	typedef Matrix<double, Dynamic, Dynamic> Mnd;
	typedef Matrix<double, 6, Dynamic> Mad;
	typedef Matrix<double, Dynamic, 6> Mad_;
	typedef Matrix<double, 6, Dynamic> Mb;
	typedef DiagonalMatrix<double, 6>  M6dd;
	typedef Matrix<double, Dynamic, 1> Vnd;
	typedef Matrix<double, 3, 1> V3d;
	typedef Matrix<double, 2, 1> V2d;
	typedef Matrix<int, Dynamic, 1> Vxi;

	typedef VectorXd Vxd;
	typedef MatrixXd Mxd;

	int IntegrationSteps;
	int PreIntegrationSteps;
	const char* ShapeType;

	double Rho;
	double L;
	double E;
	double Nu;
	double Mu;
	double A11;

	Shapes Phi;

	double Vg, Ve;
	V3d Hamiltonian;
	Vnd gradHq, gradHp, dissR;

	V7d g;
	V6d eta, deta;
	Vnd q, dq, ddq;
	Mnd M, C, Mt;
	Mad Jb, Jbt;
	Mad_ Me, Ce;
	Vnd G;

	Mb Ba;
	Mnd Sa;
	M6d Mtt, Ktt, Dtt;
	Mnd Mee, Kee, Dee;
	V6d gvec, Xi0;

	Cosserat();

	void setup(const char* str);

	void build(
		Vnd v, 	
		Vnd dv);

	void buildLagrangianModel(
		Vnd v, 
		Vnd dv);

	void buildInertiaMatrix(
		Vnd v);

	void buildGlobalMatrices();		

	void globalSystemMatODE(
		double s, 
		Mnd &K, 
		Mnd &D);

	void lagrangianODE(
		double s, 
		V19d x, 
		Mad J, 
		Mad Jt,
		V19d &dx, 
		Mad &dJ, 
		Mad &dJt, 
		Mnd &dM, 
		Mnd &dC, 
		Vnd &dG,
		Mnd &dMt, 
		Mad_ &dMe, 
		Mad_ &dCe,
		double &dVg);

	void inertiaODE(
		double s, 
		V7d x, 
		Mad J, 
		V7d &dx, 
		Mad &dJ, 
		Mnd &dM);
};

#endif

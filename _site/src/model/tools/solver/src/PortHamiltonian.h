#ifndef PORT_HAMILTONIAN_H
#define PORT_HAMILTONIAN_H

#include <unistd.h>
#include <Eigen/Dense>
#include "PortController.cpp"
#include "Cosserat.cpp"
#include "readFile.h"
#include "CARE.cpp"
//#include "forward.hpp"
//#include "eigen.hpp"

using namespace std;
using namespace Eigen;

class PortHamiltonian 
{
	public:

	int Na, Nm, Ny, Nu;

	typedef Matrix<double, Dynamic, Dynamic> Mnd;
	typedef Matrix<double, Dynamic, 1> Vnd;

	bool EC1;   // Energy controller boolean

	double h;   // time-stepping
	double t;	// time instance
	double T;	// simulation time
	Vnd q; 		// generalized coordinates
	Vnd p; 		// generalized momenta

	Vnd q_;		// coordinates estimates
	Vnd p_;		// momenta estimates
	Vnd z;		// measurement
	Vnd z_;		// measurement vel

	double H;	// Hamiltonian
	Vnd dHdq;   // Hamiltonian-Grad q
	Vnd dHdp;   // Hamiltonian-Grad p

	Mnd A; 		// linearized A SS-matrix
	Mnd B; 		// linearized B SS-matrix
	Mnd C; 		// linearized C SS-matrix

	Mnd K;		// kalman gain

	Cosserat rod;
	PortController control;

	PortHamiltonian(
		const char* str);


	void SolverUpdate();
	void ObserverUpdate();
	void TimeUpdate();

	void StromerVerletUpdateStep(
		Vnd Q0,
		Vnd P0,
		Vnd &Q1,
		Vnd &P1);

	void ImplicitTrapzoidStep(
		Vnd Q0,
		Vnd P0,
		Vnd &Q1,
		Vnd &P1);

	void LinearizedStep(
		Vnd Q0,
		Vnd P0,
		Vnd &Q1,
		Vnd &P1);

	void HessianInverse(
		Vnd R, 
		Vnd &dr);

	void ComputeJacobian(
		Vnd Q0,
		Vnd P0);

};

#endif
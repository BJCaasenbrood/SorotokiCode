#ifndef PORT_HAMILTONIAN_H
#define PORT_HAMILTONIAN_H

#include <Eigen/Dense>
#include "PortController.cpp"
#include "Cosserat.cpp"
#include "readFile.h"
//#include "forward.hpp"
//#include "eigen.hpp"

using namespace std;
using namespace Eigen;

class PortHamiltonian 
{
	public:

	int Na, Nm;

	typedef Matrix<double, Dynamic, Dynamic> Mnd;
	typedef Vector<double, Dynamic> Vnd;

	bool EC1;   // Energy controller boolean

	double h;   // time-stepping
	double t;	// time instance
	double T;	// simulation time
	Vnd q; 		// generalized coordinates
	Vnd p; 		// generalized momenta

	double H;	// Hamiltonian
	Vnd dHdq;   // Hamiltonian-Grad q
	Vnd dHdp;   // Hamiltonian-Grad p

	Cosserat rod;
	PortController control;

	PortHamiltonian(
		const char* str);

	void SolverUpdate();

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

	void HessianInverse(
		Vnd R, 
		Vnd &dr);

};

#endif
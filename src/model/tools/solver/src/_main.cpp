//https://joelcfd.com/automatic-differentiation/
#define TICTOC
#include <iostream>
#include "PortHamiltonian.cpp"
#include "TicToc.h"
#include "Logger.h"

#define DEBUG

ofstream file;
ofstream fileH;
ofstream fileG;
ofstream fileXd;
ofstream fileEta;

IOFormat matlab(5, 0, ", ", "\n", "[", "]");

void logsetup(){
	file.open("log/state.txt", ofstream::out | ofstream::trunc);
	file.close();
	fileH.open("log/hamiltonian.txt", ofstream::out | ofstream::trunc);
	fileH.close();
	fileG.open("log/endeffector.txt", ofstream::out | ofstream::trunc);
	fileG.close();
	fileXd.open("log/setpoint.txt", ofstream::out | ofstream::trunc);
	fileXd.close();
	fileEta.open("log/endeffector_Vel.txt", ofstream::out | ofstream::trunc);
	fileEta.close();
}


void showInfo(PortHamiltonian &sys){

	cout << "M(q=0): inertia matrix" << endl;
	cout << sys.rod.M.format(matlab) << endl;
	cout << endl;

	cout << "K(q=0): stiffness matrix" << endl;
	cout << sys.rod.Kee.format(matlab) << endl;
	cout << endl;

	cout << "D(q=0): stiffness matrix" << endl;
	cout << sys.rod.Dee.format(matlab) << endl;
	cout << endl;

}

int main(int argc, char** argv)
{

	// assigning rod model to PH-system/controller	
  	PortHamiltonian ph(argv[1]);

  	// clean and re-open log file
  	logsetup();
	file.open("log/state.txt", ios_base::app | ios::binary);
	fileH.open("log/hamiltonian.txt", ios_base::app | ios::binary);
	fileG.open("log/endeffector.txt", ios_base::app | ios::binary);
	fileXd.open("log/setpoint.txt", ios_base::app | ios::binary);
	fileEta.open("log/endeffector_Vel.txt", ios_base::app | ios::binary);
	//filedEta.open("log/endeffector_Acc.txt", ios_base::app | ios::binary);

	#ifdef TICTOC
		tic();
  	#endif

	#ifdef TICTOC
		showInfo(ph);
	#endif 

	// solve time integration
	while (ph.t < ph.T){

		// advance time instance
		ph.SolverUpdate();

		// output intermediate data
		fwrite(file,ph.t,ph.q);
		fwrite(fileH,ph.t,ph.rod.Hamiltonian);
		fwrite(fileG,ph.t,ph.rod.g);
		fwrite(fileXd,ph.t,ph.control.gd);
		fwrite(fileEta,ph.t,ph.rod.eta);

	}

  	#ifdef TICTOC
		toc(ph.T);
		cout <<"Frequency rate: "<< 1.0/ph.h << " Hz\n";
  	#endif

  	return EXIT_SUCCESS;
}
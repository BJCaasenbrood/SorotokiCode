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

void logsetup(){
	file.open("log/state.txt", ofstream::out | ofstream::trunc);
	file.close();
	fileH.open("log/hamiltonian.txt", ofstream::out | ofstream::trunc);
	fileH.close();
	fileG.open("log/endeffector.txt", ofstream::out | ofstream::trunc);
	fileG.close();
	fileXd.open("log/setpoint.txt", ofstream::out | ofstream::trunc);
	fileXd.close();
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

	#ifdef TICTOC
		tic();
  	#endif

	cout << ph.rod.Sa << endl;

	// solve time integration
	while (ph.t < ph.T){

		// advance time instance
		ph.SolverUpdate();

		// output intermediate data
		fwrite(file,ph.t,ph.q);
		fwrite(fileH,ph.t,ph.rod.Hamiltonian);
		fwrite(fileG,ph.t,ph.rod.g);
		fwrite(fileXd,ph.t,ph.control.gd);

	}

  	#ifdef TICTOC
		toc(ph.T);
		cout <<"Frequency rate: "<< 1.0/ph.h << " Hz\n";
  	#endif

  	return EXIT_SUCCESS;
}

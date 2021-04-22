//https://joelcfd.com/automatic-differentiation/
#define TICTOC
#include <iostream>
#include "PortHamiltonian.cpp"
#include "TicToc.h"
#include "Logger.h"

#define DEBUG

ofstream file;
ofstream file_;
ofstream fileH;
ofstream fileG;
ofstream fileXd;
ofstream fileEta;
ofstream filedEta;

IOFormat matlab(5, 0, ", ", "\n", "[", "]");

void logsetup(){
	file.open("log/state.txt", ofstream::out | ofstream::trunc);
	file.close();
	file_.open("log/estimate.txt", ofstream::out | ofstream::trunc);
	file_.close();
	fileH.open("log/hamiltonian.txt", ofstream::out | ofstream::trunc);
	fileH.close();
	fileG.open("log/endeffector.txt", ofstream::out | ofstream::trunc);
	fileG.close();
	fileXd.open("log/setpoint.txt", ofstream::out | ofstream::trunc);
	fileXd.close();
	fileEta.open("log/endeffector_Vel.txt", ofstream::out | ofstream::trunc);
	fileEta.close();
	filedEta.open("log/endeffector_Acc.txt",ofstream::out | ofstream::trunc);
	filedEta.close();
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

  	int n = ph.Na*ph.Nm;

  	// clean and re-open log file
  	logsetup();
	file.open("log/state.txt", ios_base::app | ios::binary);
	file_.open("log/estimate.txt", ios_base::app | ios::binary);
	fileH.open("log/hamiltonian.txt", ios_base::app | ios::binary);
	fileG.open("log/endeffector.txt", ios_base::app | ios::binary);
	fileXd.open("log/setpoint.txt", ios_base::app | ios::binary);
	fileEta.open("log/endeffector_Vel.txt", ios_base::app | ios::binary);
	filedEta.open("log/endeffector_Acc.txt", ios_base::app | ios::binary);

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
		fwrite(file,ph.t,ph.rod.q);
		fwrite(fileH,ph.t,ph.rod.Hamiltonian);
		fwrite(fileG,ph.t,ph.rod.g);
		fwrite(fileXd,ph.t,ph.control.gd);
		fwrite(fileEta,ph.t,ph.rod.eta);
		//fwrite(filedEta,ph.t,ph.rod.deta);
		fwrite(filedEta,ph.t,ph.control.H);
		
		ph.ObserverUpdate();
		fwrite(file_,ph.t,ph.q_);

		ph.TimeUpdate();
	}

	//cout << ph.C.block(0,n,n,1)*ph.rod.gradHp << endl;

  	#ifdef TICTOC
		toc(ph.T);
		cout <<"Frequency rate: "<< 1.0/ph.h << " Hz\n";
  	#endif

  	return EXIT_SUCCESS;
}
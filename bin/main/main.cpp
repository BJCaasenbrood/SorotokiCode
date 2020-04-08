//#include "src/RTdebug.h"
#include <fstream>
#include "src/tictoc.h"
#include "src/Model.cpp"

int main()
{
	// setting
	V6i tab;
	Vxf x;
	tab << 0,1,1,0,0,0;

	// generate model-class
	Model mdl(tab,3);

	mdl.P1 << 0,0,0.75;
	mdl.P2 << .5,0,0;

	// solve soft robotic system
	tic();
	x =  mdl.implicit_solve();
	toc();

	// output data
	ofstream myfile;
  	myfile.open ("data.log");
  	myfile << x << "\n";
  	myfile.close();

	return 0;
}


// REFERENCES
// [1] - https://github.com/stulp/tutorials/blob/master/test.md
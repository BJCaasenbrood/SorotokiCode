#include "src/Model.cpp"
//#include "src/shapesx.cpp"

int main(int argc, char** argv)
{
	// Shapes Phi;

	// Vxf stab(8), tab(2);

	// tab << 1,0;

	// stab.noalias() = tab.replicate(1,4);


	// cout << stab << endl;


	// Mxf P(2,4*2);

	// Phi.set(4,2,2,"chebyshev");

	// Phi.eval(atof(argv[1]),P);

	// cout << atof(argv[1]) << endl;

	// cout << P << endl;

	//generate model-class
	Model mdl(argv[1]);

	// // solve system
	mdl.implicit_simulate();
		
	return 0;
}

// REFERENCES
// [1] - https://github.com/stulp/tutorials/blob/master/test.md
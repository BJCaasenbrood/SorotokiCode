#include "src/Model.cpp"

int main()
{
	// setting
	V6i tab;
	tab << 0,1,1,0,0,0;

	// generate model-class
	Model mdl(tab);

	// solve system
	mdl.implicit_simulate();
	
	return 0;
}

// REFERENCES
// [1] - https://github.com/stulp/tutorials/blob/master/test.md
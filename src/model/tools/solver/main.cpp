#include "src/Model.cpp"

int main()
{
	// generate model-class
	Model mdl("config.txt");

	// solve system
	mdl.implicit_simulate();
	
	return 0;
}

// REFERENCES
// [1] - https://github.com/stulp/tutorials/blob/master/test.md
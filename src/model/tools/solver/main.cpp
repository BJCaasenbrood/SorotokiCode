#define EIGEN_NO_DEBUG
#define IMPLICIT
//#define DISCONTINIOUS
#define SOLVER_OUTPUT
#define TICTOC

#include "src/Model_ph.cpp"
using namespace std;

int main(int argc, char** argv)
{
	//generate model-class
	Model mdl(argv[1]);

	if(mdl.KINEMATIC_CONTROLLER){
		mdl.inverse_kinematics();
	}
	else if(mdl.ENERGY_CONTROLLER){
		#ifdef IMPLICIT
			mdl.implicit_simulate();
		#else
			mdl.simulate();
		#endif
	}
	else{
		#ifdef IMPLICIT
			mdl.implicit_simulate();
		#else
			mdl.simulate();
		#endif
	}

	return 0;
}

// REFERENCES
// [1] - https://github.com/stulp/tutorials/blob/master/test.md
// [2] - https://studywolf.wordpress.com/category/robotics/dynamic-movement-primitive/
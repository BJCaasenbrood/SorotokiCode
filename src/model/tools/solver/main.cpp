//#define EIGEN_NO_DEBUG
//#define DISCONTINIOUS
#define SOLVER_OUTPUT
#define TICTOC

#include "src/Model_ph.cpp"
using namespace std;

V7f g0;

int main(int argc, char** argv)
{
	//generate model-class
	Model mdl(argv[1]);

	if(mdl.KINEMATIC_CONTROLLER){
		mdl.inverse_kinematics();
		cout << "qd=" << mdl.qd.transpose().format(matlab) << endl;
	}
	else if(mdl.ENERGY_CONTROLLER){
		mdl.implicit_simulate();
	}
	else{
		//mdl.tau = mdl.q;
		//mdl.q.setZero();
		mdl.implicit_simulate();
	}

	return 0;
}

// REFERENCES
// [1] - https://github.com/stulp/tutorials/blob/master/test.md
// [2] - https://studywolf.wordpress.com/category/robotics/dynamic-movement-primitive/
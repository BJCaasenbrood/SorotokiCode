#ifndef SHAPES_H
#define SHAPES_H

#include "Eigen/Dense"

using namespace std;
using namespace Eigen;
typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;

class Shapes{

	bool cheby, poly, cubic;
  	int  NMode, NDof;

  public:
  	void set(int nmode = 3, int ndof = 6, string str = "polynomial");
  	Mxf eval(float X = 0.0);
  	Vxf phi(float X = 0.0);  	
  	
};

#endif
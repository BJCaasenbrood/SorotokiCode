#ifndef SHAPES_H
#define SHAPES_H

#include "Eigen/Dense"

using namespace std;
using namespace Eigen;
typedef Eigen::VectorXf Vxf;
typedef Eigen::MatrixXf Mxf;

class Shapes{

	bool cheby, poly, cubic, legen;
  	int  NMode, NDof;
  	float Normalize;

  public:

  	Shapes();

  	void set(int nmode = 3, int ndof = 6, const char* str = "chebyshev");
  	void eval(float X, Mxf &Phi);
  	void setNorm(float x);
  	void phi(float X, Vxf &p);  	
};

#endif
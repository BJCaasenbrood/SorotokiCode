#ifndef SHAPES_H
#define SHAPES_H

#include <iostream>
#include <Eigen/Dense>
#include "LieAlgebra.h"
#include "DEBUGGER.h"

using namespace std;
//using namespace Eigen;

class Shapes{

	bool cheby, poly, cubic, legen;
  	int  NMode, NDof;
  	double Normalize;

  public:

  	Shapes();

  	void set(int nmode = 3, int ndof = 6, const char* str = "legendre");
  	void eval(double X, Eigen::MatrixXd &Phi);
  	void setNorm(double x);
  	void phi(double X, Eigen::VectorXd &p);  	
};

#endif
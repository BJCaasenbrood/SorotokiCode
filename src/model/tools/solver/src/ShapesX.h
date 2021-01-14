#ifndef SHAPES_H
#define SHAPES_H

#include <iostream>
#include <Eigen/Dense>
#include "LieAlgebra.h"
#include "DEBUGGER.h"

using namespace std;

class Shapes{

  bool cheby, poly, cubic, legen;
  int  NMode, NDof, Nx;
  double Normalize;

  public:

    Shapes();

  	void set(int nmode = 3, int ndof = 6, 
  			 int nx = 2, const char* str = "legendre");
  	void setNorm(double x);
  	void eval(double X, Eigen::MatrixXd &Phi);
  	void phi(double X, Eigen::VectorXd &p);  	
  	
};

#endif
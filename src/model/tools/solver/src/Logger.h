#ifndef VECTOR_LOGGER_H
#define VECTOR_LOGGER_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <unistd.h>

#include <Eigen/Dense>
using namespace Eigen;

void fwrite(ofstream &file, double t, VectorXd x){

	  file << t;

	  for (int i = 0; i < x.size(); ++i)
	  	file << ", " << x(i);

	  file << "\n";
}


#endif
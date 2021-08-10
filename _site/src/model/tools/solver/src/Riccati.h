#ifndef RICCATI_H
#define RICCATI_H

#include <unistd.h>
#include <algorithm>
#include <cmath>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
using namespace Eigen;

typedef Eigen::VectorXd Vxd;
typedef Eigen::MatrixXd Mxd;

void SchurComplement(Mxd A, Mxd B) {

	int N = A.rows();
	Mxd H(2*N,2*N);

	H.block(0,0,N,N).noalias() = A;
	H.block(N,N,N,N).noalias() = -A.transpose();
	H.block(0,N,N,N).noalias() = -B*B.transpose();
	H.block(N,0,N,N).noalias() = -Mxd::Identity(N,N);

	//sleep(1200);
	RealSchur<Mxd> schur(H);
	Mxd U = schur.matrixU();
	//MatrixXd T = schur.matrixT();

	int m = U.rows();
	int n = U.cols();

	Mxd U11(m/2,n/2);
	Mxd U21(m/2,n/2);
	Mxd X(m/2,n/2);

	U11.noalias() = U.block(0,0,m/2,n/2);
	U21.noalias() = U.block(m/2,0,m/2,n/2);

	U11.transposeInPlace();
	U21.transposeInPlace();

	X.noalias() = U11.partialPivLu().solve(U21);
	X.transposeInPlace();

	cout << "A = " << A << endl;
	cout << "X = " << X << endl;
	cout << "X = " << H << endl;
	cout << "C = " << B.transpose() << endl;
 	cout << A.transpose()*X + X*A  - X*B*B.transpose()*X + Mxd::Identity(N,N) << endl;

	sleep(1200);
}

#ifndef PROPERTY_TENSOR_H
#define PROPERTY_TENSOR_H

//---------------------------------------------------
//--------------------------- build stiffness tensor
//---------------------------------------------------
void buildStiffnessTensor(
	double E,
	double Nu,
	double J11,
	double J22,
	double J33,
	double A,
	Eigen::Matrix<double, 6, 6> &Ktt
)
{
	double E0 = ((double) E);
	double G0 = ((double) E)/(2*(1+((double) Nu)));

	Ktt.setZero();
	Eigen::Matrix<double, 6, 1> v;
	v << G0*J11,E0*J22,E0*J33,E0*A,G0*A,G0*A;
	Ktt.diagonal() = v;
}

//---------------------------------------------------
//----------------------------- build inertia tensor
//---------------------------------------------------
void buildInertiaTensor(
	double Rho,
	double J11,
	double J22,
	double J33,
	double A,
	Eigen::Matrix<double, 6, 6> &Mtt
)
{
	Mtt.setZero();
	Eigen::Matrix<double, 6, 1>  v;
	v << J11,J22,J33,A,A,A;
	Mtt.diagonal() = ((double)Rho)*v;
}

//---------------------------------------------------
//----------------------------- build damping tensor
//---------------------------------------------------
void buildDampingTensor(
	double Mu,
	Eigen::Matrix<double, 6, 6> Ktt,
	Eigen::Matrix<double, 6, 6> &Dtt)
{
	Dtt = Mu*Ktt;
}

#endif
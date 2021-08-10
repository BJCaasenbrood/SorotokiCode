#include "Shapes.h"
#include <cmath>

using namespace std;
using namespace Eigen;

//---------------------------------------------------
//--------------------------------- initialize class
//---------------------------------------------------
Shapes::Shapes(){
	poly = false;
	cheby = false;
	cubic = false;
	legen = false;
}

//---------------------------------------------------
//--------------------------------- initialize class
//---------------------------------------------------
void Shapes::setNorm(double x){
	Normalize = x;
}


//---------------------------------------------------
//--------------------------------- initialize class
//---------------------------------------------------
void Shapes::set(int nmode, int ndof, const char* str) 
{
	NMode = max(0,min(nmode,10));
	NDof = max(0,min(ndof,6));

	if (str == "polynomial"){
		poly = true;
	}
	if (str == "chebyshev"){
		cheby = true;
	}
	if (str == "cubic"){
		cubic = true;
	}
	if (str == "legendre"){
		legen = true;
	}

	debug("shapes - assigned NMode and NDof");
}

//---------------------------------------------------
//---------------------------------- evaluate shapes
//---------------------------------------------------
void Shapes::phi(double s, Eigen::VectorXd &p)
{
	p.setZero();

	// construct basic polynomials
	if (poly){
		for (int i = 0; i < NMode; i++)
		{
			p(i) = pow(s,(double)(i));
		}
	}

	// construct chebyshev polynomials
	if (cheby){

		s = (2.0*s-1.0);
		for (int i = 0; i < NMode; i++)
		{
			p(i) = cos((double)(i)*acos(s));
		}
	}

		// construct legendre polynomials
	if (legen){

		s = (2.0*s-1.0);
		for (int i = 0; i < NMode; i++)
		{	
			if(i==0){ p(i) = 1;};
			if(i==1){ p(i) = s;};
			if(i==2){ p(i) = 0.5*(3*s*s - 1);};
			if(i==3){ p(i) = 0.5*(5*s*s*s - 3*s);};
			if(i==4){ p(i) = 0.125*(35*s*s*s*s - 30*s*s + 3);};
			if(i==5){ p(i) = 0.125*(63*s*s*s*s*s - 70*s*s*s + 15*s);};
			if(i==6){ p(i) = 0.0625*(231*s*s*s*s*s*s - 315*s*s*s*s + 105*s*s -5);};
			if(i==7){ p(i) = 0.0625*(429*s*s*s*s*s*s*s - 693*s*s*s*s*s + 315*s*s*s - 35*s);};
			if(i==8){ p(i) = 0.0078125*(6435*s*s*s*s*s*s*s*s - 12012*s*s*s*s*s*s + 6930*s*s*s*s - 1260*s*s + 35);};
		}
	}

}

//---------------------------------------------------
//---------------------------- evaluate shape-matrix
//---------------------------------------------------
void Shapes::eval(double s0, Eigen::MatrixXd &Phi)
{

	Eigen::VectorXd p(NMode);
	Phi.setZero();

	// normalize
	double s = s0/Normalize;

	// evaluate shape at sigma
	phi(s,p);

	// build shape function matrix
	for (int i = 0; i < NDof; i++)
	{	
		Phi.block(i,i*NMode,1,NMode).noalias() = p.transpose();
	}
}

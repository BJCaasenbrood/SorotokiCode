#include "shapes.h"
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
}

//---------------------------------------------------
//---------------------------------- evaluate shapes
//---------------------------------------------------
void Shapes::phi(float s, Vxf &p)
{
	p.setZero();

	// construct basic polynomials
	if (poly){
		for (int i = 0; i < NMode; i++)
		{
			p(i) = pow(s,(float)(i));
		}
	}

	// construct chebyshev polynomials
	if (cheby){

		s = (2.0*s-1.0);
		for (int i = 0; i < NMode; i++)
		{
			p(i) = cos((float)(i)*acos(s));
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
			//p(i) = legendre(i,s);
		}
	}

}

//---------------------------------------------------
//---------------------------- evaluate shape-matrix
//---------------------------------------------------
void Shapes::eval(float s, Mxf &Phi)
{

	Vxf p(NMode);
	Phi.setZero();

	// evaluate shape at sigma
	phi(s,p);

	// build shape function matrix
	for (int i = 0; i < NDof; i++)
	{	
		Phi.block(i,i*NMode,1,NMode).noalias() = p.transpose();
	}
}

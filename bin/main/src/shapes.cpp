#include "shapes.h"
#include <cmath>

using namespace std;
using namespace Eigen;

//---------------------------------------------------
//--------------------------------- initialize class
//---------------------------------------------------
void Shapes::set(int nmode, int ndof, string str) 
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

}

//---------------------------------------------------
//---------------------------------- evaluate shapes
//---------------------------------------------------
Vxf Shapes::phi(float s)
{
	Vxf p(NMode);
	p = Vxf::Zero(NMode);

	// construct basic polynomials
	if (poly){
		for (int i = 0; i < NMode; i++)
		{
			p(i) = pow(s,i);
		}
	}

	// construct chebyshev polynomials
	if (cheby){
		for (int i = 0; i < NMode; i++)
		{
			//p(i) = cos((float)(i)*acos((2.0*s)-1.0));
			p(i) = pow(s,(float)(i));
		}
	}

	return p;
}

//---------------------------------------------------
//---------------------------- evaluate shape-matrix
//---------------------------------------------------
Mxf Shapes::eval(float s)
{

	Vxf p(NMode);
	Mxf Phi(NDof,NDof*NMode);
	Phi = Mxf::Zero(NDof,NMode*NDof);

	// evaluate shape at sigma
	p = phi(s);

	// build shape function matrix
	for (int i = 0; i < NDof; i++)
	{	
		Phi.block(i,i*NMode,1,NMode) = p.transpose();
	}

	return Phi;
}

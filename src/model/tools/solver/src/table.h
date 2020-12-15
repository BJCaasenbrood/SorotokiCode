#ifndef TABLE_H
#define TABLE_H

#include "Eigen/Dense"
using namespace Eigen;

typedef Eigen::Array<int, Dynamic, 1> Vxi;

//---------------------------------------------------
//----------- convert table to active/constraint set
//---------------------------------------------------
Mxf tableConstraints(Vxi table, bool set){
	int k = 0;
	int na, N;

	N = table.rows();

	Mxf Id;
	Id = Mxf::Identity(N,N);

	if(set == true){na = (table>0).count();}
	else{na = (table==0).count();}

	Mxf B(N,na);

	// construct matrix of active DOF's
	if(set == true){
	for (int i = 0; i < N; ++i)
	{
		if (table(i) == true) 
		{
			B.block(0,k,N,1) = Id.block(0,i,N,1);
			k++;
		}
	}
	}

	// construct matrix of constraint DOF's
	if(set == false){
	for (int i = 0; i < N; ++i)
	{
		if (table(i) == false) 
		{
			B.block(0,k,N,1) = Id.block(0,i,N,1);
			k++;
		}
	}	
	}

	return B;

}

#endif
#include "Model.h"

using namespace std;
using namespace Eigen;

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
Model::Model(V6i tab, int nmode){
	Ba = tableConstraints(tab);
	NDof = Ba.cols();
	NMode = nmode;
	NState = NDof * NMode;

	Phi.set(NMode,NDof,"polynomial");

	buildInertiaTensor();
	buildStiffnessTensor();
	buildDampingTensor();
	buildGlobalSystem();

	q = Vxf::Zero(NState);
	dq = Vxf::Zero(NState);
	ddq = Vxf::Zero(NState);

	Xi0 << 0,0,0,1,0,0;
	P1  << 0,0,0;
	P2  << 0,0,0;

	cout << setprecision(PRECISION) << endl;
}	

//---------------------------------------------------
//----------------------- solve quasi-static problem
//---------------------------------------------------
Vxf Model::solve(){
	int n = NState;
	int i = 0;
	float t = 0; 
	float dt = (1.0*(TDOMAIN))/(1.0*(TIMESTEP));
	Vxf K1(n),K2(n);
	Vxf x(n),dx(n);

	x = Vxf::Zero(n);

	// do time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){
		kinematicODE(t,x,K1);
 		kinematicODE(t+(2.0/3.0)*dt, x+(2.0/3.0)*dt*K1, K2);
  		dx = (dt/4.0)*(K1+3.0*K2);

  		if(dx.norm() < TOLERANCE){
  			cout << "==== solver converged ===" << endl;
  			cout << "=========================" << endl;
  			cout << " iterations: " << i << endl;
  			cout << " tolerance:  " << dx.norm() << endl;
  			cout << "=========================" << endl;
  			i = MAX_ITER; 
  			x = x + dx;
  		}

  		t = t + dt;
  		x = x + dx;

	}

	// return solutions q* = q(t_eq)
	return x;
}

//----------------------- solve quasi-static problem
//---------------------------------------------------
Vxf Model::implicit_solve(){
	int n = NState;
	int i = 0;
	float t = 0; 
	float e0 = 0;
	float e1 = 0;
	float dt = (1.0*(TDOMAIN))/(1.0*(TIMESTEP));
	Vxf K1(n),K2(n);
	Vxf w(n),x(n),dx(n),R(n);
	Mxf F(n,n);

	x = Vxf::Zero(n);

	// compute time-dependent Hess
	Hess = Dee.inverse();

	// compute time-invariant Hessian
	F = -(Mxf::Identity(n,n) + dt*Dee.inverse()*Kee).inverse();

	// solve implicit time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){
		kinematicODE(t,x,K1);
		kinematicODE(t+dt,x,K2);
		R = 0.5*dt*(K1 + K2);
		dx = - F*R;

		while(R.norm() < TOLERANCE){
			kinematicODE(t+dt,x + dx,K2);
			R = dx - 0.5*dt*(K1 + K2);
			dx = dx - F*R;
		}

  		if(dx.norm() < TOLERANCE){
  			cout << "==== solver converged ===" << endl;
  			cout << "=========================" << endl;
  			cout << " iterations: " << i << endl;
  			cout << " tolerance:  " << dx.norm() << endl;
  			cout << "=========================" << endl;
  			i = MAX_ITER; 
  		}

  			t = t + dt;
  			x = x + dx;

  		i++;
	}

	// return solutions q* = q(t_eq)
	return x;
}

//---------------------------------------------------
//--------------------------- build Jacobian matrix
//---------------------------------------------------
Mxf Model::buildJacobian(float se){

	int n = NState;
	float s = 0; 
	float ds = (se)/(1.0*(SPACESTEP));

	Mxf J(6,NState);
	Mxf K1J(6,NState), K2J(6,NState);
	V7f K1,K2;
	V7f x,dx;

	x << 1,0,0,0,0,0,0;
	J = Mxf::Zero(6,NState);

	// do time integration
	for (int i = 0; i < SPACESTEP; i++){
		jacobiODE(s,x,K1,K1J);
 		jacobiODE(s+(2.0/3.0)*ds, x+(2.0/3.0)*ds*K1, K2, K2J);

  		s = s + ds;
  		x = x + (ds/4.0)*(K1+3.0*K2);
  		J = J + (ds/4.0)*(K1J+3.0*K2J);
	}

	// transform Jacobian to local frame
	return (Admap(x)).transpose()*J;
}

//---------------------------------------------------
//----------------------------------- kinematics ode
//---------------------------------------------------
void Model::kinematicODE(float t, Vxf x, Vxf &dx){
	
	int n = NMode;
	Vxf Qa, Qu;
	Mxf J1, J2, H;

	// solve static forces
	Qa = inverseDynamics(x,x*0,x*0);

	// compute manipulator Jacobian
	J1 = buildJacobian(0.5);
	J2 = buildJacobian(1.0);

	// compute pressure forces
	H = pressureMapping();
	Qu = (J1.transpose()*H)*P1;
	Qu = Qu + ((J2 - J1).transpose()*H)*P2;

	dx = Hess*(-Qa - Qu - Kee*x);
}

//---------------------------------------------------
//---------------------- inverse dynamics algorithm
//---------------------------------------------------
Vxf Model::inverseDynamics(Vxf v, Vxf dv, Vxf ddv){
	int n = 25 + NState;
	Vff x;
	Vff K1f,K2f;
	double ds,s;

	// set global states
	q = v;
	dq = dv;
	ddq = ddv;

	// init forward states
	x = Vff::Zero(19);
	x(0) = 1.0;

	// compute forward integration step
	ds = (1.0*(SDOMAIN))/(1.0*(SPACESTEP));
	s = 0.0;

	// do forward integration
	for (int i = 0; i < SPACESTEP; i++){
		forwardODE(s,x,K1f);
 		forwardODE(s+(2.0/3.0)*ds, x+(2.0/3.0)*ds*K1f, K2f);
  		x = x + (ds/4.0)*(K1f+3.0*K2f);
  		s = s + ds;
	}

	Vxf K1b(n),K2b(n);
	Vxf y(n);

	// init backward states
	y = Vxf::Zero(n);
	y.block(0,0,19,1) = x.block(0,0,19,1);

	// recompute backward integration step
	ds = (-1.0)*(1.0*(SDOMAIN))/(1.0*(SPACESTEP));
	s = 1.0*(SDOMAIN);

	// do backward integration
	for (int i = 0; i < SPACESTEP; i++){
		backwardODE(s,y,K1b);
 		backwardODE(s+(2.0/3.0)*ds, y+(2.0/3.0)*ds*K1b, K2b);
  		y = y + (ds/4.0)*(K1b+3.0*K2b);
  		s = s + ds;
	}

	// return force-vector
	return y.block(25,0,NState,1);
}

//---------------------------------------------------
//------------------------------ forward integration
//---------------------------------------------------
void Model::forwardODE(float s, Vff x, Vff &dx){

	Mxf PMat;
	M6f adxi, addxi;
	V6f xi, dxi, ddxi;
	V6f eta, deta;

	V4f quat;
	V3f xiK, xiE;
	M3f R; 
	M4f A;

	// evaluate strain-field -> sigma
	PMat = Phi.eval(s);

	xi   = (Ba*PMat)*q + Xi0;
	dxi  = (Ba*PMat)*dq;
	ddxi = (Ba*PMat)*ddq;

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	eta  = x.block(7,0,6,1);
	deta = x.block(13,0,6,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

	R = quat2rot(quat);
	A = strainMapping(R*xiK);

	// precompute adjoint actions
	adxi  = admap(xi);
	addxi = admap(dxi);

	dx.block(0,0,4,1)  = (1.0/(2*quat.norm()))*A*quat;
	dx.block(4,0,3,1)  = R*xiE;
	dx.block(7,0,6,1)  = -adxi*eta + dxi;
	dx.block(13,0,6,1) = -adxi*deta - addxi*eta + ddxi;
}

//---------------------------------------------------
//----------------------------- backward integration
//---------------------------------------------------
void Model::backwardODE(float s, Vxf x, Vxf &dx){

	int n = NState;

	Mxf PMat;
	M6f adxi, addxi, adeta;
	V6f Fb;

	V6f xi,dxi,ddxi,eta,deta,lam;
	V4f quat;
	V3f xiK, xiE;
	M3f R; 
	M4f A;

	Fb = V6f::Zero(6);

	// evaluate strain-field
	PMat = Phi.eval(s);
	xi   = (Ba*PMat)*q + Xi0;
	dxi  = (Ba*PMat)*dq;
	ddxi = (Ba*PMat)*ddq;

	R = quat2rot(quat);
	A = strainMapping(R*xiK);

	// precompute adjoint actions
	adxi  = admap(xi);
	addxi = admap(dxi);
	adeta = admap(eta);

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	eta  = x.block(7,0,6,1);
	deta = x.block(13,0,6,1);
	lam  = x.block(19,0,6,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

	dx.block(0,0,4,1)  = (1.0/(2*quat.norm()))*A*quat;
	dx.block(4,0,3,1)  = R*xiE;
	dx.block(7,0,6,1)  = -adxi*eta + dxi;
	dx.block(13,0,6,1) = -adxi*deta - addxi*eta + ddxi;
	dx.block(19,0,6,1) = adxi.transpose()*lam - (adeta.transpose()*Mtt)*eta + Mtt*deta - Fb;
	dx.block(25,0,n,1) = -(Ba*PMat).transpose()*lam;
}

//---------------------------------------------------
//-------------- forward integrate compute Jacobian
//---------------------------------------------------
void Model::jacobiODE(float s, V7f x, V7f &dx, Mxf &dJ){

	V6f xi;
	V3f xiK,xiE;
	V4f quat;
	M3f R; 
	M4f A;

	// evaluate strain-field -> sigma
	xi = (Ba*Phi.eval(s))*q + Xi0;

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

	R = quat2rot(quat);
	A = strainMapping(R*xiK);

	dx.block(0,0,4,1) = (1.0/(2*quat.norm()))*A*quat;
	dx.block(4,0,3,1) = R*xiE;

	dJ = Admap(x)*(Ba*Phi.eval(s));
}

//---------------------------------------------------
//----------- convert table to active/constraint set
//---------------------------------------------------
Mxf Model::tableConstraints(V6i table){
	int k = 0;
	int na = (table>0).count();
	Mxf B(6,na);
	M6f Id;
	Id = M6f::Identity();

	// construct matrix of active DOF's
	for (int i = 0; i < 6; ++i)
	{
		if (table(i) == true) 
		{
			B.block(0,k,6,1) = Id.block(0,i,6,1);
			k++;
		}
	}

	return B;
}

//---------------------------------------------------
//--------------------------- build stiffness tensor
//---------------------------------------------------
void Model::buildInertiaTensor(){

	Mtt = M6f::Zero(6,6);

	float A  = PI*pow(RADIUS,2);
	float J1 = 0.5*PI*pow(RADIUS,4);
	float J2 = 0.25*PI*pow(RADIUS,4);
	float J3 = 0.25*PI*pow(RADIUS,4);
	
	V6f v;
	v << J1,J2,J3,A,A,A;
	Mtt.diagonal() = ((float)RHO)*v;
}

//---------------------------------------------------
//--------------------------- build stiffness tensor
//---------------------------------------------------
void Model::buildStiffnessTensor( ){

	Ktt = M6f::Zero(6,6);

	float A  = PI*pow(RADIUS,2);
	float J1 = 0.5*PI*pow(RADIUS,4);
	float J2 = 0.25*PI*pow(RADIUS,4);
	float J3 = 0.25*PI*pow(RADIUS,4);
	
	float E0 = ((float) EMOD);
	float G0 = ((float) EMOD)/(2*(1+((float) NU)));

	V6f v;
	v << G0*J1,E0*J2,E0*J3,E0*A,G0*A,G0*A;
	Ktt.diagonal() = v;
}

//---------------------------------------------------
//----------------------------- build damping tensor
//---------------------------------------------------
void Model::buildDampingTensor(){
	Dtt = (1.0*MU)*Ktt;
}

//---------------------------------------------------
//----------------------------- build global system
//---------------------------------------------------
void Model::buildGlobalSystem(){

	Mxf P,P0,P1;
	float h, s, ds;
	s = 0.0;
	ds = (1.0*(SDOMAIN))/(1.0*(INTSTEP));

	if (INTSTEP % 3 > 0){
		cout << "WARNING! INTSTEP must be a multiple of 3! "<< endl;
	}

	P0 = Phi.eval(0.0);
	P1 = Phi.eval(1.0);

	Kee = (P0.transpose()*Ba.transpose()*Ktt*Ba*P0);
	Mee = (P0.transpose()*Ba.transpose()*Mtt*Ba*P0);
	Dee = (P0.transpose()*Ba.transpose()*Dtt*Ba*P0);
	Kee += (P1.transpose()*Ba.transpose()*Ktt*Ba*P1);
	Mee += (P1.transpose()*Ba.transpose()*Mtt*Ba*P1);
	Dee += (P1.transpose()*Ba.transpose()*Dtt*Ba*P1);

	for (int i = 1; i < INTSTEP-1; i++)
	{
		s = i*ds;
		P = Phi.eval(s);
		if (i % 3 == 0){
			Kee += 2.0*(P.transpose()*Ba.transpose()*Ktt*Ba*P);
			Mee += 2.0*(P.transpose()*Ba.transpose()*Mtt*Ba*P);
			Dee += 2.0*(P.transpose()*Ba.transpose()*Dtt*Ba*P);
		}
		else{
			Kee += 3.0*(P.transpose()*Ba.transpose()*Ktt*Ba*P);
			Mee += 3.0*(P.transpose()*Ba.transpose()*Mtt*Ba*P);
			Dee += 3.0*(P.transpose()*Ba.transpose()*Dtt*Ba*P);
		}	
	}

	Kee = (3.0*ds/8.0)*Kee;
	Mee = (3.0*ds/8.0)*Mee;
	Dee = (3.0*ds/8.0)*Dee;
}

//---------------------------------------------------
//--------------- curvature-twist to matrix operator
//---------------------------------------------------
M4f Model::strainMapping(V3f k){

	M4f K;
	float k1 = k(0);
	float k2 = k(1);
	float k3 = k(2);

	K <<  0, -k1, -k2, -k3, 
	     k1,   0, -k3,  k2, 
         k2,  k3,   0, -k1, 
         k3, -k2,  k1,   0;

	return K;
}

//---------------------------------------------------
//--------------- curvature-twist to matrix operator
//---------------------------------------------------
Mxf Model::pressureMapping(){

	Mxf K(6,3);

	K <<  0,            0,            0,
	      0, -0.5*sqrt(3), -0.5*sqrt(3),
	     -1,          0.5,          0.5,
	      0,            0,            0,
	      0,            0,            0,
	      0,            0,            0;


	return (1.0*PRS_AREA)*K;
}
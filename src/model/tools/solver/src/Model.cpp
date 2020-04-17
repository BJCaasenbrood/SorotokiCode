#include "Model.h"

using namespace std;
using namespace Eigen;

ofstream statelog;
ofstream taulog;

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
Model::Model(V6i tab){
	Ba = tableConstraints(tab);
	NDof = Ba.cols();

	NState = NDof * NMODE;

	Phi.set(NMODE,NDof,"chebyshev");

	buildInertiaTensor();
	buildStiffnessTensor();
	buildDampingTensor();
	buildGlobalSystem();

	q = Vxf::Constant(NState,ATOL);
	dq = Vxf::Zero(NState);
	ddq = Vxf::Zero(NState);
	qd = Vxf::Zero(NState);
	tau = Vxf::Zero(6);

	Xi0  << 0,0,0,1,0,0;
	gvec << 0,0,0,(float)GRAVITY,0,0;

	// read input
	read("state.log",q);
	read("input.log",tau);

	#ifdef FULL_CONTROLLER
	read("state.log",qd);
	q = Vxf::Constant(NState,ATOL);
	#endif

	#ifdef CONSTRAINED_CONTROLLER
	read("state.log",qd);
	q = Vxf::Constant(NState,ATOL);
	#endif

	// clean up 
	cleanup();


	statelog.open("state.log", ios_base::app | ios::binary);

	cout << setprecision(PRECISION) << endl;
}	

//---------------------------------------------------
//-------------------------------------- output file
//---------------------------------------------------
void Model::output(const char* str, float t, Vxf x){

  	statelog << t << " ";

  	for (int i = 0; i < x.size(); i++){
  		statelog << x(i) << " ";
  	}
  	statelog << "\n";
}

//---------------------------------------------------
//---------------------------------------- read file
//---------------------------------------------------
void Model::read(const char* str, Vxf &x){

  char data[100];
  ifstream infile;
  infile.open(str);

  for (int i = 0; i < x.size(); ++i)
  {
  	infile >> data;
  	x(i) = (float)atof(data);
  }
}

//---------------------------------------------------
//--------------------------------------- clean file
//---------------------------------------------------
void Model::cleanup(){
	ofstream myfile1, myfile2;
	myfile1.open("state.log", ofstream::out | ofstream::trunc);
	myfile2.open("tau.log", ofstream::out | ofstream::trunc);
	myfile1.close();
	myfile2.close();
}

//---------------------------------------------------
//----------------------------- real-time controller
//---------------------------------------------------
void Model::realtimeController(float t, Mxf &A1, Mxf &A2, 
	Vxf Qdes, Vxf &X){

	int n = NState;
	//Mxd G(6+n,6+n), CE(n,6+n), CI(6+n,6+n);
	//Vxd x(6+n), g0(6+n), ce0(n), ci0(6+n), Gd(6+n);
	
	Mxd G(n+6,n+6), CE(n,n+6), CI(n+6,n+6), A(n,n+6);
	Vxd x(n+6), g0(n+6), ce0(n+6), ci0(n+6), Gd(6+n), b(n);

 	Mxf H(6,3);

 	//CE.setZero();
 	G.setIdentity();
 	CI.setZero();
 	g0.setZero();
 	ci0.setZero();
 	//ce0.setZero();
 	//

 	Gd.block(0,0,6,1) = Vxd::Constant(6,1.0);
 	Gd.block(6,0,n,1) = Vxd::Constant(n,1e10);
 	G.diagonal() = Gd;
 	// ce0.noalias() = Qdes.cast<double>();

  	H.noalias() = pressureMapping();

  	A.block(0,0,n,3).noalias() = (A1.transpose()*H).cast<double>();
  	A.block(0,3,n,3).noalias() = ((A2 - A1).transpose()*H).cast<double>();
  	A.block(0,6,n,n).noalias() = (Mxf::Identity(n,n)).cast<double>();
  	//A = ((A.array().abs() > 1e-3*A.norm()).select(A.array(),0.0)).cast<double>();
  	b.noalias() = Qdes.cast<double>();
  	

  	//ce0.noalias() = -Qdes.cast<double>();
  	//G += A.transpose()*A;
  	//g0 = -A.transpose()*b;

  	//A.setZero();
  	//b.setZero();

  	// solve quadratic programming
  	solve_quadprog(G, g0, A.transpose(), -b, CI, ci0, x);
  	
  	//EigenQP::quadprog(G,g0,A,b,x);
  	// write pressure output
  	// if(isnan(x.norm()) == 0){
  	X = x.block(0,0,6,1).cast<float>();
  	// }
  	//CompleteOrthogonalDecomposition<MatrixXd> cqr(A);

  	//tau = (A.transpose() * A).ldlt().solve(A.transpose() * b).cast<float>();
  	//X.setZero();
  	// = (cqr.pseudoInverse()*b).cast<float>();

  	// cout << A << endl;
  	// cout << cqr.pseudoInverse() << endl;
  	// cout << A.cast<float>()*tau << endl;
  	// cout << b << endl;
}

//---------------------------------------------------
//------------- explicit solver quasi-static problem
//---------------------------------------------------
Vxf Model::solve(){
	int n = NState;
	int i = 0;
	float t = 0.0; 
	float dt = (1.0*(TIMESTEP));
	Vxf K1(n),K2(n);
	Vxf x(n),dx(n);

	x = q;

	// do time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){

		#ifdef WRITE_OUTPUT
		output("state.log",t,x);
		output("tau.log",t,tau);
		#endif

		kinematicODE(t,x,K1);
 		kinematicODE(t+(2.0/3.0)*dt, x+(2.0/3.0)*dt*K1, K2);
  		dx = (dt/4.0)*(K1+3.0*K2);

  		if(abs(dx.norm()) < ATOL){
  			i = MAX_ITER; 
  		}

  		x.noalias() += dx;
  		t += dt;
  		i++;
	}

	// return solutions q* = q(t_eq)
	return x;
}

//---------------------------------------------------
//------------------ explicit solver dynamic problem
//---------------------------------------------------
Vxf Model::simulate(){
	int n = NState;
	int i = 0;
	float t = 0; 
	float dt = (1.0*(TIMESTEP));

	Vxf K1(2*n),K2(2*n), K3(2*n), K4(2*n);
	Vxf x(2*n),dx(2*n);

	x = Vxf::Zero(2*n);

	// do time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){

		#ifdef WRITE_OUTPUT
		output("data.log",t,x.block(0,0,n,1));
		#endif

		dynamicODE(t,x,K1);
 		dynamicODE(t+(1.0/2.0)*dt, x+(1.0/2.0)*dt*K1, K2);
 		dynamicODE(t+(1.0/2.0)*dt, x+(1.0/2.0)*dt*K2, K3);
 		dynamicODE(t+dt, x+dt*K3, K4);
  		dx = (dt/6.0)*(K1+2.0*K2+2.0*K3+K4);

  		if(isnan(dx(0))){
  			break;
  		}

  		x.noalias() += dx;
  		t += dt;
  		i++;
	}

	// return solutions q* = q(t_eq)
	return x.block(0,0,n,1);
}

//--------------------------------------------------
//------------- implicit solve quasi-static problem
//--------------------------------------------------
Vxf Model::implicit_solve(){
	int n = NState;
	int i = 0;
	float t = 0; 
	float dt = (1.0*(TIMESTEP));

	Vxf K1(n),K2(n);
	Vxf x(n),dx(n), R(n);
	Mxf F(n,n);

	x = q;

	// compute time-invariant Hessian
	F.noalias() = (Mxf::Identity(n,n) - 0.5*dt*Dee.inverse()*Kee).inverse();
	F *= (float)SPEEDUP;

	#ifdef TICTOC
	tic();
	#endif

	// solve implicit time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){

		#ifdef WRITE_OUTPUT
		output("state.log",t,x);
		output("tau.log",t,tau);
		#endif

		kinematicODE(t,x,K1);
		kinematicODE(t+dt,x,K2);
		R.noalias() = -0.5*dt*(K1+K2);
		dx.noalias() = -F*R;

		// while(abs(R.norm()) > RTOL){
		// 	kinematicODE(t+dt,x + dx,K2);
		// 	R.noalias() = dx - 0.5*dt*(K1 + K2);
		// 	dx.noalias() -= F*R;
		// }

  		if(abs(dx.norm()) < ATOL){
  			i = MAX_ITER; 
  		}

  		x.noalias() += dx;
  		t += dt;
  		i++;
	}

	#ifdef TICTOC
	toc();
	#endif


	// return solutions q* = q(t_eq)
	return x;
}

//--------------------------------------------------
//------------------ implicit solve dynamic problem
//--------------------------------------------------
Vxf Model::implicit_simulate(){
	int n = NState;
	int i = 0;
	int j = 0;
	float t = 0; 
	float dt = (1.0*(TIMESTEP));

	Vxf K1(2*n),K2(2*n);
	Vxf x(2*n),dx(2*n),R(2*n);
	Vxf dr(2*n);

	x.block(0,0,n,1) = q;
	x.block(n,0,n,1) = dq;

	#ifdef TICTOC
	tic();
	#endif

	// solve implicit time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){

		#ifdef WRITE_OUTPUT
		output("state.log",t,x.block(0,0,n,1));
		//output("tau.log",t,tau);
		#endif

		dynamicODE(t,x,K1);

		R.noalias() = -dt*K1;
		hessianInverse(dt,R,dr);
		dx.noalias() = -dr;

		#ifndef QUASINETWON
		while(abs(R.norm()) > RTOL && j < 10){
			dynamicODE(t+dt,x+dx,K2);
			R.noalias() = dx - 0.5*dt*(K1 + K2);
			hessianInverse(dt,R,dr);
			dx.noalias() -= dr;
			j++;
		}
		#endif

  		if(isnan(dx.norm())){
  			i = MAX_ITER; 
  		}

  		x.noalias() += dx;
  		t += dt;
  		i++;
  		j = 0.0;
	}

	#ifdef TICTOC
	toc((float)TDOMAIN);
	#endif

	// return solutions q* = q(t_eq)
	return x.block(0,0,n,1);
}

//---------------------------------------------------
//-------------------------------------- dynamic ode
//---------------------------------------------------
void Model::dynamicODE(float t, Vxf x, Vxf &dx){
	
	int n = NState;
	Vxf x1(n), x2(n), null(n);
	Vxf Q(n), Qa(n), Qv(n), Qu(n), Qd(n);
	Mxf J1(6,n), J2(6,n), H(6,3);

	null.setZero();
	Qu.noalias() = null;

	// extract the generalized coordinates 
	x1.noalias() = x.block(0,0,n,1);
	x2.noalias() = x.block(n,0,n,1);	

	// solve static + dynamic forces
	inverseDynamics(x1,null,null,Qa);

	// solve static + dynamic forces
	inverseDynamics(x1,x2,null,Qv);
	Qv.noalias() -= Qa;

	// compute inertia matrix
	buildInertia(x1,Mee);

	#ifdef JACOBIAN
	// compute manipulator Jacobian
	buildJacobian(0.5, J1);
	buildJacobian(1.0, J2);
	#endif

	#ifdef CONSTRAINED_CONTROLLER
	float kp = 100;
	float kd = 5*2*sqrt(kp);

	// compute control law
	Qd.noalias() = Mee*(kp*(x1 - qd) + kd*x2) - Qa - Kee*x1;
	realtimeController(t, J1, J2, Qd, tau);
	#endif

	// compute pressure forces
	H.noalias() = pressureMapping();
	Qu.noalias() = (J1.transpose()*H)*tau.block(0,0,3,1);
	Qu.noalias() += ((J2 - J1).transpose()*H)*tau.block(3,0,3,1);

	#ifdef FULL_CONTROLLER
	float kp = 25;
	float kd = 2*sqrt(kp);
	Vxf xd1(n), xd2(n);
	xd1 = qd;
	xd2 = 0.0*qd;
	Qu.noalias() = Mee*(kp*(x1-xd1) + kd*x2) - Qa - Kee*x1;
	#endif

	Q.noalias() = Qa + Qv + Qu;

	(dx.block(0,0,n,1)).noalias() = x2;
	(dx.block(n,0,n,1)).noalias() = Mee.llt().solve(-Q - Kee*x1 - Dee*x2);

}

//---------------------------------------------------
//----------------------------------- kinematics ode
//---------------------------------------------------
void Model::kinematicODE(float t, Vxf x, Vxf &dx){
	
	int n = NState;
	Vxf Qa(n), null(n);
	Mxf J1(6,n), J2(6,n), H(6,3);

	null.setZero();

	// solve static forces
	inverseDynamics(x,null,null,Qa);

	#ifdef CONSTRAINED_CONTROLLER
	// compute manipulator Jacobian
	buildJacobian(0.5, J1);
	buildJacobian(1.0, J2);

	// compute pressure forces
	H = pressureMapping();
	Qa.noalias() += (J1.transpose()*H)*tau.block(0,0,3,1);
	Qa.noalias() += ((J2 - J1).transpose()*H)*tau.block(3,0,3,1);
	#endif

	#ifdef FULL_CONTROLLER
	float kp = 1e-8;
	float kd = 0.0;
	Qa.noalias() += kp*(q - qd);
	#endif

	dx.noalias() = Dee.inverse()*(-Qa - Kee*x);
}

//---------------------------------------------------
//----------------------- inverse dynamics algorithm
//---------------------------------------------------
void Model::inverseDynamics(Vxf v, Vxf dv, Vxf ddv, Vxf &Q){
	int n = 25 + NState;
	Vff x;
	Vff K1f,K2f;
	double ds, s;

	// set global states
	q.noalias() = v;
	dq.noalias() = dv;
	ddq.noalias() = ddv;

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
  		x.noalias() += (ds/4.0)*(K1f+3.0*K2f);
  		s += ds;
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
  		y.noalias() += (ds/4.0)*(K1b+3.0*K2b);
  		s += ds;
	}

	// return force-vector
	Q.noalias() = y.block(25,0,NState,1);
}

//---------------------------------------------------
//---------------------------- build Jacobian matrix
//---------------------------------------------------
void Model::buildJacobian(float se, Mxf &J){

	int n = NState;
	double ds,s;

	Mxf K1J(6,n), K2J(6,n);
	V13f K1, K2;
	V13f x, dx;

	x.setZero();
	x(0) = 1.0;

	// compute forward integration step
	ds = (1.0*(se))/(1.0*(SPACESTEP));
	s = 0.0;

	// initial Jacobian matrix
	J.setZero();

	// do spatial integration
	for (int i = 0; i < SPACESTEP; i++){
		jacobiODE(s,x,K1,K1J);
 		jacobiODE(s+(2.0/3.0)*ds, x+(2.0/3.0)*ds*K1, K2, K2J);

  		s += ds;
  		x.noalias() += (ds/4.0)*(K1+3.0*K2);
  		J.noalias() += (ds/4.0)*(K1J+3.0*K2J);
	}

	// transform Jacobian to local frame
	K1J.noalias() = AdmapInv(x.block(0,0,7,1))*J;
	J.noalias() = K1J;

	g.noalias()   = x.block(0,0,7,1);
	eta.noalias() = x.block(7,0,6,1);
}

//---------------------------------------------------
//----------------------------- build inertia matrix
//---------------------------------------------------
void Model::buildInertia(Vxf v, Mxf &M){

	int n = NState;
	double ds,s;

	Mxf J(6,n), K1J(6,n), K2J(6,n);
	Mxf K1M(6,n), K2M(6,n);
	V7f K1, K2;
	V7f x, dx;

	// initial matrices
	J.noalias() = Mxf::Zero(6,n);
	M.noalias() = Mxf::Zero(n,n);
	x.noalias() = V7f::Zero(7);
	x(0) = 1.0;

	// set states
	q.noalias() = v;

	// compute forward integration step
	ds = (1.0*(SDOMAIN))/(1.0*(SPACESTEP));
	s = 0.0;

	// do spatial integration
	for (int i = 0; i < SPACESTEP; i++){

  		inertiaODE(s,x,J,K1,K1J,K1M);
 		inertiaODE(s+(2.0/3.0)*ds,
 			x+(2.0/3.0)*ds*K1, 
 			J+(2.0/3.0)*ds*K1J,
 			K2,K2J,K2M);

 		s += ds;
  		x.noalias() += (ds/4.0)*(K1+3.0*K2);
  		J.noalias() += (ds/4.0)*(K1J+3.0*K2J);
  		M.noalias() += (ds/4.0)*(K1M+3.0*K2M);
	}
}

//---------------------------------------------------
//------------------------------ forward integration
//---------------------------------------------------
void Model::forwardODE(float s, Vff x, Vff &dx){

	Mxf PMat(NDof,NState);
	M6f adxi, addxi;

	// evaluate strain-field -> sigma
	Phi.eval(s,PMat);

	xi.noalias() = (Ba*PMat)*q + Xi0;
	dxi.noalias() = (Ba*PMat)*dq;
	ddxi.noalias() = (Ba*PMat)*ddq;

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	eta  = x.block(7,0,6,1);
	deta = x.block(13,0,6,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

	quat2rot(quat,R);
	strainMapping(R*xiK,A);

	// precompute adjoint actions
	admap(xi,adxi);
	admap(dxi,addxi);

	(dx.block(0,0,4,1)).noalias()  = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias()  = R*xiE;
	(dx.block(7,0,6,1)).noalias()  = -adxi*eta + dxi;
	(dx.block(13,0,6,1)).noalias() = -adxi*deta - addxi*eta + ddxi;
}

//---------------------------------------------------
//----------------------------- backward integration
//---------------------------------------------------
void Model::backwardODE(float s, Vxf x, Vxf &dx){

	int n = NState;

	Mxf PMat(NDof,NState);
	M6f adxi, addxi, adeta;
	V6f Fb;

	Fb = V6f::Zero(6);

	// evaluate strain-field
	Phi.eval(s,PMat);
	xi.noalias() = (Ba*PMat)*q + Xi0;
	dxi.noalias() = (Ba*PMat)*dq;
	ddxi.noalias() = (Ba*PMat)*ddq;

	quat2rot(quat,R);
	strainMapping(R*xiK,A);

	// precompute adjoint actions
	admap(xi,adxi);
	admap(dxi,addxi);
	admap(eta,adeta);

	// decomposition configuration space
	quat.noalias() = x.block(0,0,4,1);
	eta.noalias() = x.block(7,0,6,1);
	deta.noalias() = x.block(13,0,6,1);
	lam.noalias() = x.block(19,0,6,1);
	xiK.noalias() = xi.block(0,0,3,1);
	xiE.noalias() = xi.block(3,0,3,1);

	// add gravity component
	Fb.noalias() += (Admap(x.block(0,0,7,1)).transpose())*Mtt*gvec;

	(dx.block(0,0,4,1)).noalias()  = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias()  = R*xiE;
	(dx.block(7,0,6,1)).noalias()  = -adxi*eta + dxi;
	(dx.block(13,0,6,1)).noalias() = -adxi*deta - addxi*eta + ddxi;
	(dx.block(19,0,6,1)).noalias() = adxi.transpose()*lam - (adeta.transpose()*Mtt)*eta + Mtt*deta - Fb;
	(dx.block(25,0,n,1)).noalias() = -(Ba*PMat).transpose()*lam;
}

//---------------------------------------------------
//--------------- forward integrate compute Jacobian
//---------------------------------------------------
void Model::jacobiODE(float s, V13f x, V13f &dx, Mxf &dJ){

	Mxf PMat(NDof,NState);
	M6f adxi, Adg;

	// evaluate strain-field
	Phi.eval(s,PMat);
	xi.noalias()   = (Ba*PMat)*q + Xi0;
	dxi.noalias()  = (Ba*PMat)*dq;

	// decomposition configuration space
	quat.noalias() = x.block(0,0,4,1);
	eta.noalias()  = x.block(7,0,6,1);
	xiK.noalias()  = xi.block(0,0,3,1);
	xiE.noalias()  = xi.block(3,0,3,1);

	// precompute adjoint actions
	Adg.noalias()  = Admap(x.block(0,0,7,1));
	admap(xi,adxi);

	quat2rot(quat,R);
	strainMapping(R*xiK,A);

	(dx.block(0,0,4,1)).noalias() = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias() = R*xiE;
	(dx.block(7,0,6,1)).noalias() = -adxi*eta + dxi;

	dJ.noalias() = Adg*(Ba*PMat);
}

//---------------------------------------------------
//---------------- forward integrate compute Inertia
//---------------------------------------------------
void Model::inertiaODE(float s, V7f x, Mxf J, 
	V7f &dx, Mxf &dJ, Mxf &dM){

	Mxf PMat(NDof,NState);
	M6f Adg, AdgInv;

	// evaluate strain-field
	Phi.eval(s,PMat);
	xi.noalias() = (Ba*PMat)*q + Xi0;

	// decomposition configuration space
	quat.noalias() = x.block(0,0,4,1);
	xiK.noalias()  = xi.block(0,0,3,1);
	xiE.noalias()  = xi.block(3,0,3,1);

	quat2rot(quat,R);
	strainMapping(R*xiK,A);

	// precompute adjoint actions
	Adg.noalias() = Admap(x);
	AdgInv.noalias() = AdmapInv(x);

	(dx.block(0,0,4,1)).noalias() = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias() = R*xiE;

	// compute local D-Jacobian
	dJ.noalias() = Adg*(Ba*PMat);

	// compute local D-inertia matrix
	dM.noalias() = (AdgInv*J).transpose()*Mtt*(AdgInv*J);

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
//--------------------------- compute Hessian matrix
//---------------------------------------------------
void Model::hessianInverse(float dt, Vxf R, Vxf &dr){

	int n = NState;
	Mxf S(2*n,2*n);
	Mxf Mi(n,n);

	Mi.noalias() = Mee.householderQr().solve(Mxf::Identity(n,n));

	S.block(0,0,n,n).noalias() = -0.5*dt*Mxf::Identity(n,n);
	S.block(0,n,n,n).noalias() = Mxf::Zero(n,n);
	S.block(n,n,n,n).noalias() = 0.5*dt*Mi*Dee;
	S.block(n,0,n,n).noalias() = 0.5*dt*Mi*Kee;
	S.noalias() += Mxf::Identity(2*n,2*n);

	dr.noalias() = S.partialPivLu().solve(R);
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

	//Mxf P,P0,P1;
	int n = NState;
	double h, s, ds;

	Mxf K1K(n,n), K2K(n,n);
	Mxf K1M(n,n), K2M(n,n);
	Mxf K1D(n,n), K2D(n,n);

	s = 0.0;
	ds = (1.0*(SDOMAIN))/(1.0*(INTSTEP));

	// initial  matrix
	Kee = Mxf::Zero(n,n);
	Mee = Mxf::Zero(n,n);
	Dee = Mxf::Zero(n,n);

	// do spatial integration
	for (int i = 0; i < INTSTEP; i++){
		systemMatODE(s,K1K,K1M,K1D);
 		systemMatODE(s+(2.0/3.0)*ds,K2K,K2M,K2D);

  		s += ds;
  		Kee.noalias() += (ds/4.0)*(K1K+3.0*K2K);
  		Mee.noalias() += (ds/4.0)*(K1M+3.0*K2M);
  		Dee.noalias() += (ds/4.0)*(K1D+3.0*K2D);
	}

	Mee = ((Mee.array().abs() > 1e-2*Mee.norm()).select(Mee.array(),0.0));
	Kee = ((Kee.array().abs() > 1e-2*Kee.norm()).select(Kee.array(),0.0));
	Dee = ((Dee.array().abs() > 1e-2*Dee.norm()).select(Dee.array(),0.0));

	Mee = Mee.cast<float> ();
	Kee = Kee.cast<float> ();
	Dee = Dee.cast<float> ();

}

//---------------------------------------------------
//---------------------------- build system matrices
//---------------------------------------------------
void Model::systemMatODE(float s, 
	Mxf &K, Mxf &M, Mxf &D){

	Mxf PMat(NDof,NState);

	// evaluate strain-field
	Phi.eval(s,PMat);
	K.noalias() = ((Ba*PMat).transpose())*Ktt*(Ba*PMat);
	M.noalias() = ((Ba*PMat).transpose())*Mtt*(Ba*PMat);
	D.noalias() = ((Ba*PMat).transpose())*Dtt*(Ba*PMat);
}

//---------------------------------------------------
//--------------- curvature-twist to matrix operator
//---------------------------------------------------
Mxf Model::pressureMapping(){

	Mxf K(6,3);

	K <<  0,            0,           0,
	      0, -0.5*sqrt(3), 0.5*sqrt(3),
	   -1.0,          0.5,         0.5,
	      0,            0,           0,
	      0,            0,           0,
	      0,            0,           0;


	return ((float)PRS_AREA)*K;
}


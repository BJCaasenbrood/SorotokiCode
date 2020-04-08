#include "Model.h"

using namespace std;
using namespace Eigen;

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

	q = Vxf::Zero(NState);
	dq = Vxf::Zero(NState);
	ddq = Vxf::Zero(NState);

	Xi0  << 0,0,0,1,0,0;
	gvec << 0,0,0,(float)GRAVITY,0,0;

	cout << setprecision(PRECISION) << endl;
}	

//---------------------------------------------------
//-------------------------------------- output file
//---------------------------------------------------
void Model::output(const char* str, float t, Vxf x){

	ofstream myfile;
	myfile.open(str, ios_base::app);
  	myfile << t << " ";

  	for (int i = 0; i < x.size(); i++){
  		myfile << x(i) << " ";
  	}
  	myfile << "\n";
}

//---------------------------------------------------
//---------------------------------------- read file
//---------------------------------------------------
void Model::read(){

  tau = Vxf::Zero(6);

  char data[100];
  ifstream infile;
  infile.open("input.log");

  for (int i = 0; i < 6; ++i)
  {
  	infile >> data;
  	tau(i) = (float)atof(data);
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
void Model::realtimeController(float t, Mxf J1, Mxf J2){

	int n = NState;
	Mxf G(6,6), CE(n,6), CI(6,6);
	Vxf x(6), g0(6), ce0(6), ci0(6);
 	Mxf J0(6,n), H(6,3);
 	V6f F, u;
 	V3f gd;

 	float Kp = 1e-5;
 	float Kd = 0.5;

 	Kp *= smoothstep(t,0.0f,5.0f);

 	F.setZero();
 	G.setIdentity();
 	CI.setIdentity();
 	g0.setZero();
 	ci0.setZero();
 	gd << 0.5, 0.5, 0.5;

 	// compute control law
 	F.block(3,0,3,1).noalias() = Kp*(g.block(4,0,3,1) - gd);
 	ce0.noalias() = ((Admap(g)*J2).transpose()*F);//.cast<double>();

  	H.noalias() = pressureMapping();

  	CE.block(0,0,n,3).noalias() = (J1.transpose()*H);//.cast<double>();
  	CE.block(0,3,n,3).noalias() = ((J2 - J1).transpose()*H);//.cast<double>();
  	CE = ((CE.array().abs() > 1e-2*CE.norm()).select(CE.array(),0.0)).cast<float>();

  	// solve quadratic programming
  	//solve_quadprog(G, g0, CE.transpose(), ce0, CI, ci0, x);
	//fmincon(CE, ce0, tau);

	//cout << tau << endl;

	//sleep(1200);

  	// write pressure output
  	// = x.cast<float>();
}

//---------------------------------------------------
//------------- explicit solver quasi-static problem
//---------------------------------------------------
Vxf Model::solve(){
	int n = NState;
	int i = 0;
	float t = 0; 
	float dt = (1.0*(TDOMAIN))/(1.0*(TIMESTEP));
	Vxf K1(n),K2(n);
	Vxf x(n),dx(n);

	x = Vxf::Constant(n,ATOL);

	// do time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){

		#ifdef WRITE_OUTPUT
		output("state.log",t,x);
		output("tau.log",t,tau);
		#endif

		kinematicODE(t,x,K1);
 		kinematicODE(t+(2.0/3.0)*dt, x+(2.0/3.0)*dt*K1, K2);
  		dx = (dt/4.0)*(K1+3.0*K2);

  		if(dx.norm() < ATOL){
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
	float dt = (1.0*(TDOMAIN))/(1.0*(TIMESTEP));

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
	float dt = (1.0*(TDOMAIN))/(1.0*(TIMESTEP));

	Vxf K1(n),K2(n);
	Vxf x(n),dx(n), R(n);
	Mxf F(n,n);

	x = Vxf::Constant(n,ATOL);

	// compute time-invariant Hessian
	F.noalias() = -(Mxf::Identity(n,n) + dt*Dee.inverse()*Kee).inverse();
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
		R.noalias() = 0.5*dt*(K1 + K2);
		dx.noalias() = -F*R;

		while(abs(R.norm()) < RTOL){
			kinematicODE(t+dt,x + dx,K2);
			R.noalias() = dx - 0.5*dt*(K1 + K2);
			dx.noalias() -= F*R;
		}

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
	int i,j = 0;
	float t = 0; 
	float dt = (1.0*(TIMESTEP));

	Vxf K1(2*n),K2(2*n);
	Vxf x(2*n),dx(2*n),R(2*n);
	Mxf F(2*n,2*n);

	x = Vxf::Constant(2*n,ATOL);

	#ifdef TICTOC
	tic();
	#endif

	// solve implicit time integration
	while (t < ((float)TDOMAIN) && i < MAX_ITER){

		#ifdef WRITE_OUTPUT
		output("state.log",t,x.block(0,0,n,1));
		output("tau.log",t,tau);
		#endif

		dynamicODE(t,x,K1);
		dynamicODE(t+dt,x,K2);

		R.noalias() = 0.5*dt*(K1 + K2);
		F.noalias() = -hessianInverse(dt);
		dx.noalias() = -F*R;

		while(abs(R.norm()) > RTOL){
			dynamicODE(t+dt,x+dx,K2);
			R.noalias() -= 0.5*dt*(K1 + K2);
			
			#ifdef QUASINETWON
			F.noalias() = -hessianInverse(dt);
			#endif
			
			dx.noalias() -= F*R;
			j++;
		}

  		if(isnan(dx.norm())){
  			//i = MAX_ITER; 
  		}

  		//cout << j << endl;

  		x.noalias() += dx;
  		t += dt;
  		j = 0;
  		i++;
	}

	#ifdef TICTOC
	toc((float)TDOMAIN);
	#endif

	// return solutions q* = q(t_eq)
	return x.block(0,0,n,1);
}

//--------------------------------------------------
//------------ min f(x) := ||Ax - b||2 with x >= 0
//--------------------------------------------------
void Model::fmincon(Mxf &A, Vxf &b, Vxf &x){

	int n = x.size();
	int i = 0;
	float lambda = 1.0;
	float R = 1.0e3;
	Mxf F,I;

	Vxf r(n),d(n);

	// initialize
	d.setZero();
	r.setZero();

	// add numerical robustness
	A = ((A.array().abs() > 1e-2*A.norm()).select(A.array(),0.0)).cast<float>();
	I = Mxf::Identity(A.rows(),A.cols());
	F.noalias() = (lambda*A.adjoint()*A + I).inverse();

	while(abs(R) > ATOL && i < MAX_ITER){
		// optimize routine
		d.noalias() = pthres(x-r);
		x.noalias() = F*(pthres(d) + r + lambda*A.adjoint()*b);
		r.noalias() = r + pthres(x) - x;

		// check residual
		R = (A*x - b + pthres(-x)).norm();
		i++;
	}

	//cout << i << endl;

}

//---------------------------------------------------
//-------------------------------------- dynamic ode
//---------------------------------------------------
void Model::dynamicODE(float t, Vxf x, Vxf &dx){
	
	int n = NState;
	Vxf x1(n), x2(n), null(n);
	Vxf Qa(n), Qu(n), Qv(n);
	Mxf M(n,n), J1(6,n), J2(6,n), H(6,3);

	// extract the generalized coordinates 
	x1 = x.block(0,0,n,1);
	x2 = x.block(n,0,n,1);	
	null = Vxf::Zero(n);

	// solve static + dynamic forces
	Qa.noalias() = inverseDynamics(x1,x2,null);

	// compute manipulator Jacobian
	buildJacobian(0.5, J1);
	buildJacobian(1.0, J2);

	// compute control law
	//realtimeController(t, J1, J2);

	// compute pressure forces
	H.noalias() = pressureMapping();
	Qa.noalias() += (J1.transpose()*H)*tau.block(0,0,3,1);
	Qa.noalias() += ((J2 - J1).transpose()*H)*tau.block(3,0,3,1);

	// compute inertia matrix
	buildInertia(x1,Mee);

	(dx.block(0,0,n,1)).noalias() = x2;
	(dx.block(n,0,n,1)).noalias() = (Mee).inverse()*(-Qa - Kee*x1 - Dee*x2);
}

//---------------------------------------------------
//----------------------------------- kinematics ode
//---------------------------------------------------
void Model::kinematicODE(float t, Vxf x, Vxf &dx){
	
	int n = NState;
	Vxf Qa(n), Qu(n), null(n);
	Mxf J1(6,n), J2(6,n), H(6,3);

	null = Vxf::Zero(n);

	// solve static forces
	Qa.noalias() = inverseDynamics(x,null,null);

	// compute manipulator Jacobian
	buildJacobian(0.5, J1);
	buildJacobian(1.0, J2);

	// compute pressure forces
	H = pressureMapping();
	Qu.noalias() = (J1.transpose()*H)*tau.block(0,0,3,1);
	Qu.noalias() += ((J2 - J1).transpose()*H)*tau.block(3,0,3,1);

	dx.noalias() = Dee.inverse()*(-Qa - Qu - Kee*x);
}

//---------------------------------------------------
//----------------------- inverse dynamics algorithm
//---------------------------------------------------
Vxf Model::inverseDynamics(Vxf v, Vxf dv, Vxf ddv){
	int n = 25 + NState;
	Vff x;
	Vff K1f,K2f;
	double ds, s;

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
	return y.block(25,0,NState,1);
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
	J = ((Admap(x.block(0,0,7,1))).transpose())*J;

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
	J = Mxf::Zero(6,n);
	M = Mxf::Zero(n,n);
	x = V7f::Zero(7);
	x(0) = 1.0;
	q = v;

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

	Mxf PMat;
	M6f adxi, addxi;

	// evaluate strain-field -> sigma
	PMat = Phi.eval(s);

	xi.noalias() = (Ba*PMat)*q + Xi0;
	dxi.noalias() = (Ba*PMat)*dq;
	ddxi.noalias() = (Ba*PMat)*ddq;

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	eta  = x.block(7,0,6,1);
	deta = x.block(13,0,6,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

	R.noalias() = quat2rot(quat);
	A.noalias() = strainMapping(R*xiK);

	// precompute adjoint actions
	adxi.noalias() = admap(xi);
	addxi.noalias() = admap(dxi);

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

	Mxf PMat;
	M6f adxi, addxi, adeta;
	V6f Fb;

	Fb = V6f::Zero(6);

	// evaluate strain-field
	PMat.noalias() = Phi.eval(s);
	xi.noalias() = (Ba*PMat)*q + Xi0;
	dxi.noalias() = (Ba*PMat)*dq;
	ddxi.noalias() = (Ba*PMat)*ddq;

	R.noalias() = quat2rot(quat);
	A.noalias() = strainMapping(R*xiK);

	// precompute adjoint actions
	adxi.noalias() = admap(xi);
	addxi.noalias() = admap(dxi);
	adeta.noalias() = admap(eta);

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	eta  = x.block(7,0,6,1);
	deta = x.block(13,0,6,1);
	lam  = x.block(19,0,6,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

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

	Mxf PMat;
	M6f adxi, Adg;

	// evaluate strain-field
	PMat.noalias() = Phi.eval(s);
	xi.noalias()   = (Ba*PMat)*q + Xi0;
	dxi.noalias()  = (Ba*PMat)*dq;

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	eta  = x.block(7,0,6,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

	// precompute adjoint actions
	Adg.noalias()  = Admap(x.block(0,0,7,1));
	adxi.noalias() = admap(xi);

	R.noalias() = quat2rot(quat);
	A.noalias() = strainMapping(R*xiK);

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

	Mxf PMat;
	M6f Adg;

	// evaluate strain-field
	PMat.noalias() = Phi.eval(s);
	xi.noalias() = (Ba*PMat)*q + Xi0;

	// decomposition configuration space
	quat = x.block(0,0,4,1);
	xiK  = xi.block(0,0,3,1);
	xiE  = xi.block(3,0,3,1);

	R.noalias() = quat2rot(quat);
	A.noalias() = strainMapping(R*xiK);

	// precompute adjoint actions
	Adg.noalias() = Admap(x);

	(dx.block(0,0,4,1)).noalias() = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias() = R*xiE;

	dJ.noalias() = Adg*(Ba*PMat);

	dM.noalias() = ((J.transpose())*Adg)*Mtt*((Adg.transpose())*J);
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
Mxf Model::hessianInverse(float dt){

	int n = NState;
	Mxf F(2*n,2*n);
	F.block(0,0,n,n).noalias() = dt*Mxf::Identity(n,n);
	F.block(0,n,n,n).noalias() = dt*Mxf::Zero(n,n);
	F.block(n,n,n,n).noalias() = dt*Mee.inverse()*Dee;
	F.block(n,0,n,n).noalias() = dt*Mee.inverse()*Kee;
	F.noalias() += Mxf::Identity(2*n,2*n);

	return F.inverse();
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

	Mxf P;

	// evaluate strain-field
	P = Phi.eval(s);
	K.noalias() = ((Ba*P).transpose())*Ktt*(Ba*P);
	M.noalias() = ((Ba*P).transpose())*Mtt*(Ba*P);
	D.noalias() = ((Ba*P).transpose())*Dtt*(Ba*P);
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

	K <<  0,            0,           0,
	      0, -0.5*sqrt(3), 0.5*sqrt(3),
	     -1,          0.5,         0.5,
	      0,            0,           0,
	      0,            0,           0,
	      0,            0,           0;


	return ((float)PRS_AREA)*K;
}


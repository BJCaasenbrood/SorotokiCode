#include "Model.h"

using namespace std;
using namespace Eigen;

ofstream statelog;	// log of states
ofstream taulog;	// log of control inputs
ofstream glog;		// log of configuration
ofstream elog;		// log of energy variables

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
Model::Model(const char* str){

	V6i tab;

	ConfigFile cf(str);

	ENERGY_CONTROLLER = static_cast<bool>(cf.Value("options","ENERGY_CONTROLLER"));
	WRITE_OUTPUT = static_cast<bool>(cf.Value("options","WRITE_OUTPUT"));

	K1 = static_cast<int>(cf.Value("cosserat","K1"));	// e_xy
	K2 = static_cast<int>(cf.Value("cosserat","K2"));	// e_xz
	K3 = static_cast<int>(cf.Value("cosserat","K3"));	// e_yz
	E1 = static_cast<int>(cf.Value("cosserat","E1"));	// e_xx
	E2 = static_cast<int>(cf.Value("cosserat","E2"));	// e_yy
	E3 = static_cast<int>(cf.Value("cosserat","E3"));	// e_zz

	tab << K1,K2,K3,E1,E2,E3;

	NDISC   = static_cast<int>(cf.Value("model","NDISC"));
	NMODE   = static_cast<int>(cf.Value("model","NMODE"));
  	SDOMAIN = cf.Value("model","SDOMAIN");
  	TDOMAIN = cf.Value("model","TDOMAIN");

  	SPACESTEP = static_cast<int>(cf.Value("solver","SPACESTEP"));
  	INTSTEP   = static_cast<int>(cf.Value("solver","INTSTEP"));
  	MAX_IMPL  = static_cast<int>(cf.Value("solver","MAX_IMPL"));
  	MAX_ITER  = static_cast<int>(cf.Value("solver","MAX_ITER"));
  	ATOL      = cf.Value("solver","ATOL");
  	RTOL      = cf.Value("solver","RTOL");
  	SPEEDUP   = cf.Value("solver","SPEEDUP");
  	TIMESTEP  = cf.Value("solver","TIMESTEP");

  	RHO      = cf.Value("physics","RHO");
  	EMOD     = cf.Value("physics","EMOD");
  	NU       = cf.Value("physics","NU");
  	MU       = cf.Value("physics","MU");
  	PRS_AREA = cf.Value("physics","PRS_AREA");
  	GRAVITY  = cf.Value("physics","GRAVITY");
  	RADIUS   = cf.Value("physics","RADIUS");

  	KP = cf.Value("control","KP");
  	KD = cf.Value("control","KD");
  	KE = cf.Value("control","KE");

	Ba = tableConstraints(tab,true);
	Bc = tableConstraints(tab,false);
	NDof = Ba.cols();
	NState = NDof * NMODE;

	Vxi sa(NDISC), stab(NState);
	sa.setZero();
	sa(0) = 1;
	stab = sa.replicate(NMODE,1);

	Sa = tableConstraints(stab,true);
	Sc = tableConstraints(stab,false);

	Sa.transposeInPlace();
	Sc.transposeInPlace();

	Phi.set(NMODE,NDof,NDISC,"chebyshev");

	buildInertiaTensor();
	buildStiffnessTensor();
	buildDampingTensor();
	buildGlobalSystem();

	q   = Vxf::Constant(NState,ATOL);
	dq  = Vxf::Zero(NState);
	ddq = Vxf::Zero(NState);
	qd  = Vxf::Zero(NState);
	tau = Vxf::Zero(NState);
	u   = Vxf::Zero(NMODE);
	z0  = Vxf::Zero(6);

	Xi0 << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
	gvec << 0.0, 0.0, 0.0, GRAVITY, 0.0, 0.0;

	// read input
	read("state.log",q);
	read("state_dt.log",dq);
	read("input.log",tau);
	read("point.log",z0);

	#ifdef FULL_CONTROLLER
		read("state.log",qd);
		qd = Vxf::Constant(NState,ATOL);
	#endif

	#ifdef CONSTRAINED_CONTROLLER
		read("state.log",qd);
		q = Vxf::Constant(NState,ATOL);
	#endif

	// if(ENERGY_CONTROLLER){
	// 	read("state.log",qd);

	// 	// compute desired energy variable;
	// 	inverseDynamics(qd, qd*0, qd*0, q);
	// 	Ec = Hamil(0) + Hamil(1);

	// 	cout << "Desired H:=T+V: " << Ec << endl;

	// 	q = Vxf::Constant(NState,ATOL);
	// }

	// clean up 
	cleanup();

	statelog.open("state.log", ios_base::app | ios::binary);
	taulog.open("tau.log", ios_base::app | ios::binary);
	glog.open("g.log", ios_base::app | ios::binary);
	elog.open("e.log", ios_base::app | ios::binary);

	cout << setprecision(PRECISION) << endl;
}	

//---------------------------------------------------
//-------------------------------------- output file
//---------------------------------------------------
void Model::output(ofstream &file, float t, Vxf x){

	  file << t << " ";

	  for (int i = 0; i < x.size(); ++i)
	  {
	  	file << x(i) << " ";
	  }

	  file << "\n";
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
	ofstream myfile1, myfile2, myfile3, myfile4;
	myfile1.open("state.log", ofstream::out | ofstream::trunc);
	myfile2.open("tau.log", ofstream::out | ofstream::trunc);
	myfile3.open("g.log", ofstream::out | ofstream::trunc);
	myfile4.open("e.log", ofstream::out | ofstream::trunc);
	myfile1.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
}

// //--------------------------------------------------
// //----------- transformation to task-space dynamics
// //--------------------------------------------------
// void Model::operationalSpaceDynamics(Mxf &J, Mxf &Jt, 
// 	Vxf &dq, Mxf &M, Vxf &C, Vxf &G, Mxf &Mx, Vxf &Cx, 
// 	Vxf &Gx){

// 	int n = NState;
// 	int k = Mx.rows();

// 	Mxf Mi(n,n), Li(k,k), Jg(n,k);

// 	// compute task-space inertia matrix
// 	Mi.noalias() = M.householderQr().solve(Mxf::Identity(n,n));
// 	Li.noalias() = J*(Mi*J.transpose());
// 	Mx.noalias() = Li.householderQr().solve(Mxf::Identity(k,k));
		
// 	// build generalized inverse
// 	Jg.noalias() = Mi*J.transpose()*Mx;
// 	//Jg.transposeInPlace();
	
// 	// residual Lagrangian forces in task-space
// 	Cx.noalias() = Jg.transpose()*C - Mx*Jt*dq;
// 	Gx.noalias() = Jg.transpose()*G;	

// }

// //--------------------------------------------------
// //----- dynamically consistent null-space projector
// //--------------------------------------------------
// void Model::dynamicProjector(Mxf &J, Mxf &M, Mxf &S,
// 	Mxf &P){

// 	int n = NState;
// 	int m = J.rows();
// 	int k = S.cols();

// 	Mxf Mi(n,n), X(m,k), Xpv(k,m);

// 	// precompute inverse inertia
// 	Mi.noalias() = M.householderQr().solve(Mxf::Identity(n,n));
// 	Xpv.noalias() = (J*Mi*S).completeOrthogonalDecomposition().pseudoInverse();

// 	// build dynamically consistent projector
// 	P.noalias() = Xpv*J*Mi;

// }
// //---------------------------------------------------
// //----------- controller wrench at end-effector level
// //---------------------------------------------------
// void Model::controllerWrench(float t, Mxf &J, Vxf &f){

// 	V6f z, dz, dz0, tmp;
// 	M6f Adg;

// 	dz0.setZero();

// 	// get configuration-space in R6
// 	SE3toR6(g,z);
// 	Admap(g,Adg);

// 	dz.noalias() = eta;
// 	z0.block(0,0,3,1).noalias() = z.block(0,0,3,1);
// 	dz0.block(0,0,3,1).noalias() = eta.block(0,0,3,1);

// 	// compute wrench
// 	f = KP*(z - z0) + KD*(dz - dz0);
// 	f *= smoothstep(0.2*t,0,1);
// }

//---------------------------------------------------
//-- passivity-based controller for end-effector
//---------------------------------------------------
void Model::controllerPassive(float t, Mxf &M, Vxf &Qv,
Vxf &Qa, Mxf &J, Vxf &f){

	int na,nc;
	na = NDof*NDISC;
	nc = NState - na;

	//V6f z, dz, dz0;
	Vxf u_(na);
	//M6f Adg;

	//dz0.setZero();

	// // get configuration-space in R6
	//SE3toR6(g,z);
	//Admap(g,Adg);

	//dz.noalias() = eta;
	//z0.block(0,0,3,1).noalias() = z.block(0,0,3,1);
	//dz0.block(0,0,3,1).noalias() = eta.block(0,0,3,1);

	// // compute wrench
	//u_ = Sa*J.transpose()*(KP*(z - z0) + KD*(dz - dz0));
	// f *= smoothstep(0.1*t,0,1);

	Mxf M11(na,na), M12(na,nc);
	Mxf M21(nc,na), M22(nc,nc); 
	Mxf Ms(NState,NState);
	Mxf M22_(na,na);
	Vxf H2_(na);
	Mxf S(NState,NState);

	S.block(0,0,na,NState).noalias()  = Sc;
	S.block(na,0,nc,NState).noalias() = Sa;

	// partition of inertia matrix;
	Ms.noalias() = S*M*S.transpose();
	M11.noalias() = Ms.block(0,0,nc,nc);
	M12.noalias() = Ms.block(0,nc,nc,na);
	M21.noalias() = Ms.block(nc,0,na,nc);
	M22.noalias() = Ms.block(nc,nc,na,na);

	M22_.noalias() = M22 - M21*M11.householderQr().solve(M12);
	H2_.noalias() = Sa*Qa - M21*M11.householderQr().solve(Sc*Qa);

	// energy shaping
	float E_;
	E_ = (Hamil(0) + Hamil(1)) - Ec;

	f.noalias() = M22_*(-KP*Sa*(q - qd) - KD*Sa*dq) + H2_;
	f *= smoothstep(t,0,1);
}

/*
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

*/

/*
//---------------------------------------------------
//------------------ explicit solver dynamic problem
//---------------------------------------------------
//
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

		dynamicODE(t,x,K1);
 		dynamicODE(t+(1.0/2.0)*dt, x+(1.0/2.0)*dt*K1, K2);
 		dynamicODE(t+(1.0/2.0)*dt, x+(1.0/2.0)*dt*K2, K3);
 		dynamicODE(t+dt, x+dt*K3, K4);
  		dx = (dt/6.0)*(K1+2.0*K2+2.0*K3+K4);

  		if(isnan(dx(0))){
  			break;
  		}

  		#ifdef WRITE_OUTPUT
		output(statelog,t,x.block(0,0,n,1));
		output(taulog,t,u);
		output(glog,t,g);
		#endif

		#ifdef SOLVER_OUTPUT
		cout << "-------------/ " << endl;
		cout << "tim =" << t << endl;
		#endif 

  		x.noalias() += dx;
  		t += dt;
  		i++;
	}

	// return solutions q* = q(t_eq)
	return x.block(0,0,n,1);
}
*/

/*
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
*/

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
	while (t < TDOMAIN && i < MAX_ITER){

		dynamicODE(t,x,K1);

		R.noalias() = -dt*K1;
		hessianInverse(dt,R,dr);
		dx.noalias() = -dr;

		#ifndef QUASINETWON
		while(abs(R.block(n,0,n,1).norm()) > RTOL && j < MAX_IMPL){
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

  		if(WRITE_OUTPUT){
			output(statelog,t,x.block(0,0,n,1));
			output(taulog,t,u);
			output(glog,t,g);
			output(elog,t,Hamil);
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
	Mxf J1t(6,n), J2t(6,n);

	null.setZero();
	Qu.setZero();

	// extract the generalized coordinates 
	x1.noalias() = x.block(0,0,n,1);
	x2.noalias() = x.block(n,0,n,1);	

	// solve position-dependent forces
	inverseDynamics(x1,null,null,Qa);

	// solve velocity-dependent forces
	inverseDynamics(x1,x2,null,Qv);
	Qv.noalias() -= Qa;

	// add elastic material forces
	Qa.noalias() += Kee*x1;

	// add viscous material forces
	Qv.noalias() += Dee*x2;

	// compute inertia matrix
	buildInertia(x1,Mee);

	// compute pressure forces
	H.noalias() = pressureMapping();

	// compute manipulator Jacobian
	#ifdef JACOBIAN
		buildJacobian(0.5, J1, J1t);
		buildJacobian(1.0, J2, J2t);
	#endif

	if(ENERGY_CONTROLLER){
		int l = Ba.cols();
		int m = 4; // dimension of input-space
		int k = 6; // dimension of task-space

		//Mxf S(n,m), Prj(k,n);
		Mxf J(k,n), Jt(k,n);
		//Mxf Mz(k,k);
		//Vxf Cz(k), Pz(k), f(6);
		//Vxf u(NDof*NDISC);
		M6f Adg;

		// recover task-space Jacobian 
		Admap(g,Adg);
		J.noalias()  = Adg*J2;
		Jt.noalias() = Adg*J2t;

		//S.block(0,0,n,l).noalias() = J1.transpose()*Ba;
  		//S.block(0,l,n,l).noalias() = (J2 - J1).transpose()*Ba;
  		//S.noalias() = J2.transpose()*Ba;

  		// compute operational space dynamics
		//operationalSpaceDynamics(J,Jt,x2,Mee,Qv,Qa,Mz,Cz,Pz);

		// compute null space projection
		//dynamicProjector(J,Mee,S,Prj);

		// compute control-action in operational space
		//controllerWrench(t,J,f);
		
		controllerPassive(t,Mee,Qv,Qa,J,u);

		Qu.noalias() = Sa.transpose()*u;

		// project torques
		//Qu.noalias() = S*Prj*J.transpose()*(f);
		//u.noalias() = Prj*J.transpose()*(f);
	}

	//Qu.noalias() += (J1.transpose()*H)*tau.block(0,0,3,1);
	//Qu.noalias() += ((J2 - J1).transpose()*H)*tau.block(3,0,3,1);
	//
	//Qu.noalias() += ((float)PRS_AREA)*tau;

	//sleep(1200);

	#ifdef FULL_CONTROLLER
		float kp = 25;
		float kd = 2*sqrt(kp);
		Vxf xd1(n), xd2(n);
		xd1 = qd;
		xd2 = 0.0*qd;
		Qu.noalias() = Mee*(kp*(x1-xd1) + kd*x2) - Qa - Kee*x1;
	#endif

	(dx.block(0,0,n,1)).noalias() = x2;
	(dx.block(n,0,n,1)).noalias() = Mee.llt().solve(-Qa - Qv + Qu);

}

//---------------------------------------------------
//----------------------------------- kinematics ode
//---------------------------------------------------
void Model::kinematicODE(float t, Vxf x, Vxf &dx){
	
	int n = NState;
	Vxf Qa(n), null(n);
	Mxf J1(6,n), J2(6,n), H(6,3);
	Mxf J1t(6,n), J2t(6,n);

	null.setZero();

	// solve static forces
	inverseDynamics(x,null,null,Qa);

	#ifdef CONSTRAINED_CONTROLLER
	// compute manipulator Jacobian
	buildJacobian(0.5, J1, J1t);
	buildJacobian(1.0, J2, J2t);

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
	int n = 27 + NState;
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

	// write Hamiltonian function (minus for back-integration)
	Hamil.noalias() = -y.block(25,0,2,1);

	// return force-vector
	Q.noalias() = y.block(27,0,NState,1);
}

//---------------------------------------------------
//---------------------------- build Jacobian matrix
//---------------------------------------------------
void Model::buildJacobian(float se, Mxf &J, Mxf &Jt){

	int n = NState;
	double ds,s;

	Mxf K1J(6,n), K2J(6,n);
	Mxf K1Jt(6,n), K2Jt(6,n);
	M6f adeta, AdgInv;
	V13f K1, K2;
	V13f x, dx;

	// initial Jacobian matrix
	J.setZero();
	Jt.setZero();
	x.setZero();
	x(0) = 1.0;

	// compute forward integration step
	ds = (1.0*(se))/(1.0*(SPACESTEP));
	s = 0.0;

	// do spatial integration
	for (int i = 0; i < SPACESTEP; i++){
		jacobiODE(s,x,K1,K1J,K1Jt);
 		jacobiODE(s+(2.0/3.0)*ds, x+(2.0/3.0)*ds*K1, K2, K2J, K2Jt);

  		s += ds;
  		x.noalias() += (ds/4.0)*(K1+3.0*K2);
  		J.noalias() += (ds/4.0)*(K1J+3.0*K2J);
  		Jt.noalias() += (ds/4.0)*(K1Jt+3.0*K2Jt);
	}

	// return configuration and velocities
	g.noalias()   = x.block(0,0,7,1);
	eta.noalias() = x.block(7,0,6,1);

	// compute adjoint actions
	AdmapInv(g,AdgInv);
	admap(eta,adeta);

	// transform Jacobian to local frame
	K1J.noalias() = AdgInv*J;
	J.noalias() = K1J;

	// transform time-derivative Jacobian to local frame
	K1Jt.noalias() = AdgInv*Jt;
	Jt.noalias()  = K1Jt;
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
	M6f adxi, addxi, adeta, Adg;
	V6f Fg;

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
	Admap(x.block(0,0,7,1),Adg);

	// decomposition configuration space
	quat.noalias() = x.block(0,0,4,1);
	eta.noalias() = x.block(7,0,6,1);
	deta.noalias() = x.block(13,0,6,1);
	lam.noalias() = x.block(19,0,6,1);
	xiK.noalias() = xi.block(0,0,3,1);
	xiE.noalias() = xi.block(3,0,3,1);

	// add gravity component
	Fg.noalias() = (Adg).transpose()*Mtt*gvec;

	(dx.block(0,0,4,1)).noalias()  = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias()  = R*xiE;
	(dx.block(7,0,6,1)).noalias()  = -adxi*eta + dxi;
	(dx.block(13,0,6,1)).noalias() = -adxi*deta - addxi*eta + ddxi;
	(dx.block(19,0,6,1)).noalias() = adxi.transpose()*lam - (adeta.transpose()*Mtt)*eta + Mtt*deta - Fg;
	(dx.block(25,0,1,1)).noalias() = 0.5*(eta.transpose()*Mtt)*eta;	// kinetic energy
	(dx.block(26,0,1,1)).noalias() = lam.transpose()*xi;			// potential energy
	(dx.block(27,0,n,1)).noalias() = -(Ba*PMat).transpose()*lam;	// force projection
}

//---------------------------------------------------
//--------------- forward integrate compute Jacobian
//---------------------------------------------------
void Model::jacobiODE(float s, V13f x, V13f &dx, Mxf &dJ, Mxf &dJt){

	Mxf PMat(NDof,NState);
	M6f adxi, Adg, adeta;

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
	Admap(x.block(0,0,7,1),Adg);
	admap(xi,adxi);
	admap(eta,adeta);

	quat2rot(quat,R);
	strainMapping(R*xiK,A);

	(dx.block(0,0,4,1)).noalias() = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias() = R*xiE;
	(dx.block(7,0,6,1)).noalias() = -adxi*eta + dxi;

	dJ.noalias() = Adg*(Ba*PMat);
	dJt.noalias() = Adg*adeta*(Ba*PMat);
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
	Admap(x,Adg);
	AdmapInv(x,AdgInv);

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
Mxf Model::tableConstraints(Vxi table, bool set){
	int k = 0;
	int na,N;

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
void Model::buildStiffnessTensor(){

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

	Mee = ((Mee.array().abs() > 1e-6*Mee.norm()).select(Mee.array(),0.0));
	Kee = ((Kee.array().abs() > 1e-6*Kee.norm()).select(Kee.array(),0.0));
	Dee = ((Dee.array().abs() > 1e-6*Dee.norm()).select(Dee.array(),0.0));

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
	   -0.0,         -0.0,        -0.0,
	      0,            0,           0,
	      0,            0,           0;


	return ((float)PRS_AREA)*K;
}


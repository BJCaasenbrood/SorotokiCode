#include "Model_ph.h"

using namespace std;
using namespace Eigen;

IOFormat matlab(5, 0, ", ", "\n", "[", "]");

ofstream statelog_dt;
ofstream statelog;
ofstream taulog;
ofstream glog;
ofstream xdlog;

//--------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
Model::Model(const char* str){

	V6i tab;

	int E1,E2,E3;
  	int K1,K2,K3; 
  	float Xd, Yd, Zd;
  	float Q1d, Q2d, Q3d, Q4d;
  	int ACTUATION_SPACE;

	ConfigFile cf(str);

	WRITE_OUTPUT = static_cast<bool>(cf.Value("options","WRITE_OUTPUT"));
	ENERGY_CONTROLLER = static_cast<bool>(cf.Value("options","ENERGY_CONTROLLER"));
	KINEMATIC_CONTROLLER = static_cast<bool>(cf.Value("options","KINEMATIC_CONTROLLER"));
	ACTUATION_SPACE = static_cast<int>(cf.Value("options","ACTUATION_SPACE"));
	
	if(KINEMATIC_CONTROLLER & ENERGY_CONTROLLER){
		cout << "error: both controllers active. Switchting to kinematic controller." << endl;
		ENERGY_CONTROLLER = false;
	}

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
  	MAX_IK    = static_cast<int>(cf.Value("solver","MAX_IK"));
  	LAMBDA    = cf.Value("solver","LAMBDA");
  	ATOL      = cf.Value("solver","ATOL");
  	RTOL      = cf.Value("solver","RTOL");
  	SPEEDUP   = cf.Value("solver","SPEEDUP");
  	TIMESTEP  = cf.Value("solver","TIMESTEP");
  	ADAMP     = cf.Value("solver","ADAMPING");

  	RHO      = cf.Value("physics","RHO");
  	EMOD     = cf.Value("physics","EMOD");
  	NU       = cf.Value("physics","NU");
  	MU       = cf.Value("physics","MU");
  	PRS_AREA = cf.Value("physics","PRS_AREA");
  	GRAVITY  = cf.Value("physics","GRAVITY");
  	RADIUS   = cf.Value("physics","RADIUS");
  	A11      = cf.Value("physics","AREA");
  	J11      = cf.Value("physics","J_XX");
  	J22      = cf.Value("physics","J_YY");
  	J33      = cf.Value("physics","J_ZZ");

  	KP = cf.Value("control","KP");
  	KD = cf.Value("control","KD");
  	Xd = cf.Value("setpoint","Xd");
  	Yd = cf.Value("setpoint","Yd");
  	Zd = cf.Value("setpoint","Zd");
  	Q1d = cf.Value("setpoint","Q1d");
  	Q2d = cf.Value("setpoint","Q2d");
  	Q3d = cf.Value("setpoint","Q3d");
  	Q4d = cf.Value("setpoint","Q4d");

	Ba = tableConstraints(tab,true);
	Bc = tableConstraints(tab,false);
	NDof = Ba.cols();
	NState = NDof * NMODE;

	Vxi sa(NMODE/NDISC), sc(NMODE), stab(NState);
	stab.setZero();
	stab(0) = 1;

	// building actuation space
	#ifdef DISCONTINIOUS
		if (ACTUATION_SPACE == -1){
			sa.setConstant(1.0);
		}
		else{
			sa(0) = 1;
		}

		sc = sa.replicate(NDISC,1);
		Phi.set(NMODE,NDof,NDISC,"legendre");
		
	#else
		if (ACTUATION_SPACE == -1){
			sc.setConstant(1.0);
		}
		else{
			sc.setZero();
			sc(0) = 1.0;
		}

		Phi.set(NMODE,NDof,"legendre");
		
	#endif

	// normalize domain
	Phi.setNorm(SDOMAIN);

	stab = sc.replicate(NDof,1);
	Sa = tableConstraints(stab,true);
	Sc = tableConstraints(stab,false);

	Sa.transposeInPlace();
	Sc.transposeInPlace();

	buildInertiaTensor();
	buildStiffnessTensor();
	buildDampingTensor();
	buildGlobalSystem();

	q    = Vxf::Constant(NState,ATOL);
	dq   = Vxf::Zero(NState);
	ddq  = Vxf::Zero(NState);
	dq_  = Vxf::Zero(NState);
	qd   = Vxf::Zero(NState);
	tau  = Vxf::Zero(NState);
	u    = Vxf::Zero(6);
	z0   = Vxf::Zero(6);
	gvec = Vxf::Zero(6);

	Xi0 << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
	gd << Q1d, Q2d, Q3d, Q4d, Xd, Yd, Zd;

	// read input files
	q    = readVecXf("log/state.log");
	dq   = readVecXf("log/state_dt.log");
	gvec = readVecXf("log/grav_vector.log");
	tau  = readVecXf("log/tau_vector.log");

	// clean up 
	cleanup();

	statelog.open("log/state.log", ios_base::app | ios::binary);
	statelog_dt.open("log/statelog_dt.log", ios_base::app | ios::binary);
	taulog.open("log/tau.log", ios_base::app | ios::binary);
	glog.open("log/g.log", ios_base::app | ios::binary);
	xdlog.open("log/xd.log", ios_base::app | ios::binary);

	// pre-build lagrangian system
	buildLagrange(q,dq,Mee,Cee,Gee,Mtee,Mbee,Cbee);

	// compute desired Hamiltonian
	// dq(0) = 1.0;
	// Hstar = Hamiltonian(dq,ddq);
	// dq.setZero();

	// setting output precision
	cout << setprecision(PRECISION) << endl;
}	

//---------------------------------------------------
//-------------------------------------- output file
//---------------------------------------------------
void Model::output(ofstream &file, float t, Vxf x){

	  file << t;

	  for (int i = 0; i < x.size(); ++i)
	  file << ", " << x(i);

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
	ofstream file;

	file.open("log/state.log", ofstream::out | ofstream::trunc);
	file.close();

	file.open("log/tau.log", ofstream::out | ofstream::trunc);
	file.close();
	
	file.open("log/g.log", ofstream::out | ofstream::trunc);
	file.close();

	file.open("log/statelog_dt.log", ofstream::out | ofstream::trunc);
	file.close();

	file.open("log/xd.log", ofstream::out | ofstream::trunc);
	file.close();
}

//---------------------------------------------------
//-- solve minimal potential energy kinematics
//---------------------------------------------------
void Model::inverse_kinematics(float t){

	float dt = clamp(SPEEDUP*(1.0*(TIMESTEP)),1e-4,0.99);

	Mxf J(6,NState);
	Mxf Jt(6,NState);
	Vxf dx(NState);

	Mxf V,H;
	Cxd cod;

	Mxf N(NState,NState);
	M4f G,Gi;
	M6f T;
	V6f Phi, etad;

	// store current q in qd 
	qd.noalias() = q;
	dx(0) = 1;

	int i = 0;
	while (i < MAX_IK && abs(dx.norm()) > ATOL ){

		// compute jacobian matrix
		buildJacobian((float)SDOMAIN, J, Jt);

		// compute eta desired
		SE3(g,G);
		SE3Inv(gd,Gi);
		logmapSE3(Gi*G,Phi);
		tmapSE3(Phi,T);
		etad.noalias() = T*Phi;

		// compute kernel of J
		cod.compute(J);
		V.noalias() = cod.matrixZ().transpose();
		H.noalias() = cod.colsPermutation()*V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());

		N.noalias() = H*(H.transpose()*Dee*H).partialPivLu().solve(H.transpose());

		dx.noalias() = -N*((Sc.transpose()*Sc)*Kee*q) - (Mxf::Identity(NState,NState) - N*Dee)*J.transpose()*
		(J*J.transpose() + LAMBDA*Mxf::Identity(6,6)).partialPivLu().solve(etad);

		// update states
		if(!isnan(dx.norm())){
  			q.noalias() += dt*dx;	
  		}

  		if(WRITE_OUTPUT){
  			t = i*dt;
			output(statelog,t,q);
			output(taulog,t,tau);
			output(glog,t,g);
			output(xdlog,t,gd);
		}

		i++;
	}

	#ifdef TICTOC
		toc((float)TDOMAIN);
		cout <<"Number of solver steps: "<< i << "#\n";
	#endif

	if(WRITE_OUTPUT){
		output(statelog,t,q);
		output(xdlog,t,gd);
		output(glog,t,g);
	}

	dx.noalias() = q;
	q.noalias()  = qd;
	qd.noalias() = dx;

}

//---------------------------------------------------
//-- passivity-based controller for end-effector
//---------------------------------------------------
void Model::controllerPassive(float t, Vxf &Hq, Vxf &Hp, Vxf &f){

	int n = NState;
	int m = Sa.cols();
	Mxf Xpv(m,6);

	M6f K;
	Mxf J(6,n);
	Mxf Jt(6,n);
	Vxf dx(n);
	V6f R;
	V6f k, xd;

	Mxf V,H;
	Cxd cod;

	Mxf N(NState,NState);
	M4f G,Gi;
	M6f T;
	V6f Phi, etad;

	R.setZero();
	K.setZero();
	xd.setZero();

	// set virtual stiffness gains
	k << 1e-5,1e-5,1e-5,1,1,1;
	K.diagonal() = k;

	// end-effector jacobian
	buildJacobian((float)SDOMAIN, J, Jt);

	// relative spatial-twist between g(l) and g*
	SE3(g,G);
	SE3Inv(gd,Gi);
	logmapSE3(Gi*G,Phi);
	tmapSE3(Phi,T);
	R.noalias() = K*T*Phi*smoothstep(t,0.0,5.0);

	// compute dq - desired potential energy 
	dx.noalias() = J.transpose()*(J*J.transpose() + 
		LAMBDA*Mxf::Identity(6,6)).householderQr().solve(K*R);

	// energy-shaping + damping injection
	f.noalias() = Sa*(Hq - KP*dx - KD*Dee*Hp);
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
	float alpha;

	Vxf K1(2*n),K2(2*n);
	Vxf x(2*n),dx(2*n),R(2*n);
	Vxf x1(2*n),x2(2*n);
	Vxf dr(2*n);
	V7f gtmp;

	x1.setZero();
	x2.setZero();
	dx.setZero();

	x.block(0,0,n,1) = q;
	x.block(n,0,n,1) = dq;

	#ifdef TICTOC
		tic();
	#endif

	// backup desired configuration
	gtmp.noalias() = gd;

	// hessian update parameter	
	alpha = (2.0/3.0)*dt;			

	// solve implicit time integration
	while (t < TDOMAIN && i < MAX_ITER){

		// update trajectory
		pathSolve(t,gd);

		// dynamics ode
		dynamicODE(t,x,K1);

		R.noalias() = -dt*K1;
		hessianInverse(alpha,R,dr);
		dx.noalias() = -dr;

		while(abs(R.block(0,0,n,1).norm()) > RTOL && j < MAX_IMPL){
			dynamicODE(t+dt,x+dx,K2);
			R.noalias() = dx - 0.5*dt*(K1 + K2);
			hessianInverse(alpha,R,dr);
			dx.noalias() -= dr;
			j++;
		}

  		if(isnan(dx.norm())){
  			i = MAX_ITER; 
  		}

  		if(WRITE_OUTPUT){
			output(statelog,t,x.block(0,0,n,1));
			output(statelog_dt,t,x.block(n,0,n,1));
			output(taulog,t,tau);
			output(glog,t,g);
			output(xdlog,t,gd);
		}

		// state/time update
  		x.noalias() += dx;
  		t += dt;

  		// recover d-configuration
  		gd.noalias() = gtmp;
  		j = 0.0;
	}

	#ifdef TICTOC
		toc((float)TDOMAIN);
		cout <<"Frequency rate: "<< 1.0/dt << " Hz\n";
	#endif

	// return solutions q* = q(t_eq)
	return x.block(0,0,n,1);
}

//--------------------------------------------------
//------------------ implicit solve dynamic problem
//--------------------------------------------------
Vxf Model::simulate(){
	int n = NState;
	int i = 0;
	float t = 0; 
	float dt = (1.0*(TIMESTEP));
	float alpha;

	Vxf K1(2*n), K2(2*n), K3(2*n), K4(2*n);
	Vxf x(2*n),dx(2*n);
	V7f gtmp;

	x.setZero();
	dx.setZero();

	x.block(0,0,n,1) = q;
	x.block(n,0,n,1) = dq;

	#ifdef TICTOC
		tic();
	#endif

	// backup desired configuration
	gtmp.noalias() = gd;

	// update parameter	
	alpha = (1.0/2.0)*dt;		

	// do time integration
	while (t < TDOMAIN && i < MAX_ITER){

		// update trajectory
		pathSolve(t,gd);

		dynamicODE(t,x,K1);
 		dynamicODE(t+alpha, x+alpha*K1, K2);
 		dynamicODE(t+alpha, x+alpha*K2, K3);
 		dynamicODE(t+dt, x+dt*K3, K4);

  		dx = (dt/6.0)*(K1+2.0*K2+2.0*K3+K4);

  		if(isnan(dx.norm())){
  			i = MAX_ITER; 
  		}

  		if(WRITE_OUTPUT){
			output(statelog,t,x.block(0,0,n,1));
			output(statelog_dt,t,x.block(n,0,n,1));
			output(taulog,t,tau);
			output(glog,t,g);
			output(xdlog,t,gd);
		}

  		// state/time update
  		x.noalias() += dx;
  		t += dt;

  		// recover d-configuration
  		gd.noalias() = gtmp;
	}

	#ifdef TICTOC
		toc((float)TDOMAIN);
		cout <<"Frequency rate: "<< 1.0/dt << " Hz\n";
	#endif

	// return solutions q* = q(t_eq)
	return x.block(0,0,n,1);
}

//---------------------------------------------------
//-------------------------------------- dynamic ode
//---------------------------------------------------
void Model::dynamicODE(float t, Vxf x, Vxf &dx){
	
	int n = NState;
	Vxf x1(n), x2(n);
	Vxf dHdp, dHdq;
	Vxf Q(n), Qa(n), Qu(n);

	Qa.setZero();
	Qu.setZero();

	// extract the generalized coordinates 
	x1.noalias() = x.block(0,0,n,1); 	// q-vector
	x2.noalias() = x.block(n,0,n,1);	// p-vector

	// compute part-diff p Hamiltonian
	dHdp.noalias() = Mee.householderQr().solve(x2);	

	// compute Lagrangian model
	buildLagrange(x1,dHdp,Mee,Cee,Gee,Mtee,Mbee,Cbee);

	// compute part-diff q Hamiltonian
	dHdq.noalias() = -(Mtee - Cee)*dHdp + Gee + Kee*x1;

	// compute hamiltonian/total energy
	float H = Hamiltonian(x1,x2);

	if(ENERGY_CONTROLLER){

		//tau.noalias() = KP*(H - Hstar)*Sa*dHdp;
		controllerPassive(t,dHdq,dHdp,tau);
		
		Qu.noalias() = Sa.transpose()*tau;
	}

	else{

		Qu.noalias() = Sa.transpose()*tau;
	}

	(dx.block(0,0,n,1)).noalias() = dHdp;
	(dx.block(n,0,n,1)).noalias() = -dHdq - Dee*dHdp + Qu + Qa;
}

//---------------------------------------------------
//---------------------------- build Jacobian matrix
//---------------------------------------------------
void Model::buildJacobian(float se, Mxf &J, Mxf &Jt){

	int n = NState;
	double ds,s;

	Mxf K1J(6,n), K2J(6,n);
	Mxf K1Jt(6,n), K2Jt(6,n);
	M6f AdgInv;
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
  		x.noalias()  += (ds/4.0)*(K1+3.0*K2);
  		J.noalias()  += (ds/4.0)*(K1J+3.0*K2J);
  		Jt.noalias() += (ds/4.0)*(K1Jt+3.0*K2Jt);
	}

	// return configuration and velocities
	g.noalias()   = x.block(0,0,7,1);
	eta.noalias() = x.block(7,0,6,1);

	// compute adjoint actions
	AdmapInv(g,AdgInv);
	//admap(eta,adeta);

	// transform Jacobian to local frame
	K1J.noalias() = AdgInv*J;
	J.noalias() = K1J;

	// transform time-derivative Jacobian to local frame
	K1Jt.noalias() = AdgInv*Jt;
	Jt.noalias()  = K1Jt;
}

//---------------------------------------------------
//--------------------------  build lagrangian model
//---------------------------------------------------
void Model::buildLagrange(Vxf v, Vxf dv, 
	Mxf &M, Mxf &C, Vxf &G, Mxf &Mt, Mxf &Me, Mxf &Ce){

	int n = NState;
	double ds,s;

	Mxf J(6,n), K1J(6,n), K2J(6,n);
	Mxf Jt(6,n), K1Jt(6,n), K2Jt(6,n);
	Mxf K1M(n,n), K2M(n,n);
	Mxf K1Mt(n,n), K2Mt(n,n);
	Mxf K1C(n,n), K2C(n,n);
	Vxf K1G(n), K2G(n);
	Mxf K1Me(n,6), K2Me(n,6);
	Mxf K1Ce(n,6), K2Ce(n,6);
	V13f K1, K2;
	V13f x, dx;

	// initial matrices
	J.setZero();
	Jt.setZero();
	M.setZero();
	Mt.setZero();
	C.setZero();
	G.setZero();
	Ce.setZero();
	Me.setZero();
	x.setZero();
	x(0) = 1.0;

	// set states
	q.noalias()  = v;
	dq.noalias() = dv;

	// compute forward integration step
	ds = (1.0*(SDOMAIN))/(1.0*(SPACESTEP));
	s = 0.0;

	// do spatial integration
	for (int i = 0; i < SPACESTEP; i++){

  		lagrangianODE(s,x,J,Jt,K1,K1J,K1Jt,K1M,K1C,K1G,K1Mt,K1Me,K1Ce);

 		lagrangianODE(s+(2.0/3.0)*ds,
 				x+(2.0/3.0)*ds*K1, 
 				J+(2.0/3.0)*ds*K1J,
 				Jt+(2.0/3.0)*ds*K1Jt,
 				K2,K2J,K2Jt,K2M,K2C,K2G,K2Mt,K2Me,K2Ce);

 		s += ds;
  		x.noalias()  += (ds/4.0)*(K1+3.0*K2);  
  		J.noalias()  += (ds/4.0)*(K1J+3.0*K2J);
  		Jt.noalias() += (ds/4.0)*(K1Jt+3.0*K2Jt);
  		M.noalias()  += (ds/4.0)*(K1M+3.0*K2M);
  		C.noalias()  += (ds/4.0)*(K1C+3.0*K2C);
  		G.noalias()  += (ds/4.0)*(K1G+3.0*K2G);
  		Mt.noalias() += (ds/4.0)*(K1Mt+3.0*K2Mt);
  		Me.noalias() += (ds/4.0)*(K1Me+3.0*K2Me);
  		Ce.noalias() += (ds/4.0)*(K1Ce+3.0*K2Ce);

	}
}

//---------------------------------------------------
//--------------- forward integrate compute Jacobian
//---------------------------------------------------
void Model::jacobiODE(float s, V13f x, V13f &dx, 
	Mxf &dJ, Mxf &dJt){

	Mxf PMat(NDof,NState);
	M6f adxi, Adg, adeta;

	// evaluate strain-field
	Phi.eval(s,PMat);
	xi.noalias()  = (Ba*PMat)*q + Xi0;
	dxi.noalias() = (Ba*PMat)*dq;

	// decomposition configuration space
	g.noalias()    = x.block(0,0,7,1);
	quat.noalias() = x.block(0,0,4,1);
	eta.noalias()  = x.block(7,0,6,1);
	xiK.noalias()  = xi.block(0,0,3,1);
	xiE.noalias()  = xi.block(3,0,3,1);

	// precompute adjoint actions
	Admap(g,Adg);
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
//--------------- forward integrate Lagrangian model
//---------------------------------------------------
void Model::lagrangianODE(float s, V13f x, Mxf J, Mxf Jt,
	V13f &dx, Mxf &dJ, Mxf &dJt, Mxf &dM, Mxf &dC, Vxf &dG,
	Mxf &dMt, Mxf &dMe, Mxf &dCe){

	Mxf PMat(NDof, NState);
	M6f Adg, Adg_;
	M6f adxi, adeta;

	// evaluate strain-field
	Phi.eval(s,PMat);
	xi.noalias()   = (Ba*PMat)*q + Xi0;
	dxi.noalias()  = (Ba*PMat)*dq;

	// decomposition configuration space
	g.noalias()    = x.block(0,0,7,1);
	quat.noalias() = x.block(0,0,4,1);
	eta.noalias()  = x.block(7,0,6,1);
	xiK.noalias()  = xi.block(0,0,3,1);
	xiE.noalias()  = xi.block(3,0,3,1);

	quat2rot(quat,R);
	strainMapping(R*xiK,A);

	// precompute adjoint actions
	Admap(g,Adg);
	AdmapInv(g,Adg_);
	admap(xi,adxi);
	admap(eta,adeta);

	// compute local D-configuration
	(dx.block(0,0,4,1)).noalias() = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias() = R*xiE;
	(dx.block(7,0,6,1)).noalias() = -adxi*eta + dxi;

	// compute local D-Jacobian
	dJ.noalias()  = Adg*(Ba*PMat);
	dJt.noalias() = Adg*adeta*(Ba*PMat);

	// compute local D-inertia matrix
	dM.noalias() = (Adg_*J).transpose()*Mtt*(Adg_*J);

	// compute local C-coriolis matrix
	dC.noalias() = (Adg_*J).transpose()*((Mtt*adeta - 
	 	adeta.transpose()*Mtt)*(Adg_*J) + Mtt*(Adg_*Jt));

	// compute local G-potential vector
	dG.noalias() = (Adg_*J).transpose()*Mtt*Adg_*gvec;

	// compute local time-derivative of G-inertia matrix
	dMt.noalias() = (Adg_*Jt).transpose()*Mtt*(Adg_*J) + (Adg_*J).transpose()*Mtt*(Adg_*Jt);

	// compute local moving-base at \eta_0: Me and Ce matrix
	dMe.noalias() = (Adg_*J).transpose()*Mtt;
	dCe.noalias() = (Adg_*J).transpose()*(Mtt*adeta - adeta.transpose()*Mtt);
}

//---------------------------------------------------
//--------------------------- compute Hessian matrix
//---------------------------------------------------
void Model::hessianInverse(float alpha, Vxf R, Vxf &dr){

	int n = NState;
	Mxf S(2*n,2*n);
	Mxf Mi(n,n);

	Mi.noalias() = Mee.partialPivLu().solve(Mxf::Identity(n,n));

	S.block(0,n,n,n).noalias() = -alpha*Mi;
	S.block(0,0,n,n).noalias() = Mxf::Zero(n,n);
	S.block(n,0,n,n).noalias() = alpha*Kee;
	S.block(n,n,n,n).noalias() = alpha*Dee*Mi - alpha*(Mtee - Cee)*Mi;
	S.noalias() += Mxf::Identity(2*n,2*n);

	dr.noalias() = S.partialPivLu().solve(R);
}

//---------------------------------------------------
//---------------- build nonlinear elastic potential
//---------------------------------------------------
void Model::buildNonlinearElastic(Vxf x, Vxf &N){

	int n = NState;;
	float Ne1,Ne2,Ne3;
	float Nb1,Nb2,Nb3;

	float ke1 = 2.23e3;
	float ke2 = 1.73e3;
	float ke3 = -4.55e2;

	float kb1 = 4.23e-1;
	float kb2 = 3.99e-1;
	float kb3 = -2.29e-1;
}

//---------------------------------------------------
//----------------------------- build inertia tensor
//---------------------------------------------------
void Model::buildInertiaTensor(){

	Mtt.setZero();

	V6f v;
	v << J11,J22,J33,A11,A11,A11;
	Mtt.diagonal() = ((float)RHO)*v;
}

//---------------------------------------------------
//--------------------------- build stiffness tensor
//---------------------------------------------------
void Model::buildStiffnessTensor(){

	Ktt.setZero();

	float E0 = ((float) EMOD);
	float G0 = ((float) EMOD)/(2*(1+((float) NU)));

	V6f v;
	v << G0*J11,E0*J22,E0*J33,E0*A11,G0*A11,G0*A11;
	Ktt.diagonal() = v;
}

//---------------------------------------------------
//----------------------------- build damping tensor
//---------------------------------------------------
void Model::buildDampingTensor(){
	Dtt.setZero();
	Dtt = (1.0*MU)*Ktt + (0.0*MU)*Mtt;
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

	// initial  matrix
	Kee  = Mxf::Zero(NState,NState);
	Mee  = Mxf::Zero(NState,NState);
	Dee  = Mxf::Zero(NState,NState);
	Cee  = Mxf::Zero(NState,NState);
	Mtee = Mxf::Zero(NState,NState);
	Mbee = Mxf::Zero(NState,6);
	Cbee = Mxf::Zero(NState,6);
	Gee  = Vxf::Zero(NState);

	s = 0.0;
	ds = (1.0*(SDOMAIN))/(1.0*(INTSTEP));

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
//------------- pressure mapping for 3dof soft robot
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

//---------------------------------------------------
//------------------------------ compute Hamiltonian
//---------------------------------------------------
float Model::Hamiltonian(Vxf x1, Vxf x2){

	float H;

	H =  0.5*x2.transpose()*Mee.householderQr().solve(x2);
	H =+ 0.5*x1.transpose()*Kee*x1;

	return H;
}

//---------------------------------------------------
//------- Poisson bracket [f,g] = (f,q*g,p - f,p*g,q)
//---------------------------------------------------
void PoissonBracket(Vxf &Hq, Vxf &Hp, Mxf &A){
	//A.noalias() = -Hp*Sa.transpose();
}

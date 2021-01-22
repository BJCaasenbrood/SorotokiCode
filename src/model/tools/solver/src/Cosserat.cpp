#include "Cosserat.h"

//#define DEBUG

//--------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
Cosserat::Cosserat()
{

}

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
void Cosserat::setup(const char* str){

	ConfigFile cf(str);

	Vxi Table(6);
	int E1,E2,E3;
  	int K1,K2,K3; 
  	int ActuationSpace;
  	int AssignProperties;

  	Na = static_cast<int>(cf.Value("model","NDOF"));
  	Nm = static_cast<int>(cf.Value("model","NMODE"));
  	Nd = static_cast<int>(cf.Value("model","NDISC"));

  	K1 = static_cast<int>(cf.Value("cosserat","K1"));	// e_xy
	K2 = static_cast<int>(cf.Value("cosserat","K2"));	// e_xz
	K3 = static_cast<int>(cf.Value("cosserat","K3"));	// e_yz
	E1 = static_cast<int>(cf.Value("cosserat","E1"));	// e_xx
	E2 = static_cast<int>(cf.Value("cosserat","E2"));	// e_yy
	E3 = static_cast<int>(cf.Value("cosserat","E3"));	// e_zz

	L   = cf.Value("physics","LENGTH");
	Rho = cf.Value("physics","RHO");
  	E   = cf.Value("physics","EMOD");
  	Nu  = cf.Value("physics","NU");
  	Mu  = cf.Value("physics","MU");

  	A11 = cf.Value("physics","AREA");
  	double J11 = cf.Value("physics","J_XX");
  	double J22 = cf.Value("physics","J_YY");
  	double J33 = cf.Value("physics","J_ZZ");

  	Mee.noalias() = MatrixXd::Zero(Na*Nm,Na*Nm);
  	Kee.noalias() = MatrixXd::Zero(Na*Nm,Na*Nm);
  	Dee.noalias() = MatrixXd::Zero(Na*Nm,Na*Nm);
  	M.noalias()   = MatrixXd::Zero(Na*Nm,Na*Nm);
  	Mt.noalias()  = MatrixXd::Zero(Na*Nm,Na*Nm);
  	C.noalias()   = MatrixXd::Zero(Na*Nm,Na*Nm);
  	Jb.noalias()  = MatrixXd::Zero(6,Na*Nm);
  	Jbt.noalias() = MatrixXd::Zero(6,Na*Nm);
  	Me.noalias()  = MatrixXd::Zero(Na*Nm,6);
  	Ce.noalias()  = MatrixXd::Zero(Na*Nm,6);
  	Ba.noalias()  = MatrixXd::Zero(6,Na);
  	G.noalias()   = VectorXd::Zero(Na*Nm);
  	q.noalias()   = VectorXd::Constant(Na*Nm,1e-3);
  	dq.noalias()  = VectorXd::Zero(Na*Nm);

  	gradHp.noalias() = VectorXd::Zero(Na*Nm);
  	gradHq.noalias() = VectorXd::Zero(Na*Nm);

	Table << K1,K2,K3,E1,E2,E3;

	AssignProperties = static_cast<int>(cf.Value("options","PROPERTYSET"));
	ActuationSpace = static_cast<int>(cf.Value("options","ACTUSPACE"));
	PreIntegrationSteps = cf.Value("solver","INTSTEP");
	IntegrationSteps    = cf.Value("solver","SPACESTEP");

	Ba.noalias() = tableConstraints(Table,true);
	debug("Table constraint build");

	Vxi sa(Nm/Nd), sc(Nm), stab(Na*Nm);
	//sa.setZero();
	stab.setZero();
	stab(0) = 1;

	#ifdef DISCONTINIOUS
		if (ActuationSpace == -1){
			sa.setConstant(1.0);
		}
		else{
			sa(0) = 1;
		}

		sc = sa.replicate(Nd,1);
		Phi.set(Nm,Na,Nd,"legendre");
		Phi.setNorm(L);
		
	#else
		if (ActuationSpace == -1){
			sc.setConstant(1.0);
		}
		else{
			sc.setZero();
			sc(0) = 1.0;
		}

		Phi.set(Nm,Na,"legendre");
		Phi.setNorm(L);

	#endif

	stab = sc.replicate(Na,1);
	Sa = tableConstraints(stab,true);
	Sa.transposeInPlace();
	
	Xi0  = readVecXf("log/xi0_vector.txt");
	gvec = readVecXf("log/grav_vector.txt");

	// precompute density-properties
	if (AssignProperties == -1){
		buildStiffnessTensor(E,Nu,J11,J22,J33,A11,Ktt);
		buildInertiaTensor(Rho,J11,J22,J33,A11,Mtt);
		buildDampingTensor(Mu,Ktt,Dtt);
	}
	else{
		buildInertiaTensor(Rho,J11,J22,J33,A11,Mtt);
		Ktt.setZero();
		Ktt.diagonal() = readVecXf("log/Ktt_vector.txt");
		buildDampingTensor(Mu,Ktt,Dtt);
	}

	// prebuild global matrices
	buildGlobalMatrices();

	// prebuild Lagrangian model
	buildLagrangianModel(q,dq);

	// prebuild PortHamiltonian model
	build(q,M*dq);
	debug("Assembly of lagrangian matrices");

}

//--------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
void Cosserat::build(
	Vnd x1, 
	Vnd x2)
{
	// compute generalized inertia matrix
	buildInertiaMatrix(x1);
	debug("building generalized inertia matrix");

	// compute part-diff p Hamiltonian = dq
	gradHp.noalias() = M.householderQr().solve(x2);	

	// compute Lagrangian model
	buildLagrangianModel(x1,gradHp);
	debug("building full lagrangian model");

	// compute nonlinear stiffness
	//buildstiffness(x1,Kee);

	// compute part-diff q Hamiltonian
	gradHq.noalias() = -(Mt - C)*gradHp + G + Kee*x1;//+ Dee*gradHp;

	// compute total energy - Hamiltonian
	Hamiltonian(0) = 0.5*x2.dot(gradHp);
	Hamiltonian(1) = 0.5*x1.dot(Kee*x1);
	Hamiltonian(2) = Vg;

	debug("Cosserat building");
}

//---------------------------------------------------
//------------------------------ build gloabl system
//---------------------------------------------------
void Cosserat::buildGlobalMatrices(
)
{
	double s, ds;
	double TOL = 1e-6;

	Mnd K1K(Na*Nm,Na*Nm), K2K(Na*Nm,Na*Nm);
	Mnd K1D(Na*Nm,Na*Nm), K2D(Na*Nm,Na*Nm);

	// initial  matrix
	Kee.setZero();
	Dee.setZero();

	s = 0.0;
	ds = (1.0*(L))/(1.0*(PreIntegrationSteps));

	// do spatial integration
	for (int i = 0; i < PreIntegrationSteps; i++){
		globalSystemMatODE(s,K1K,K1D);
 		globalSystemMatODE(s+(2.0/3.0)*ds,K2K,K2D);

  		s += ds;
  		Kee.noalias() += (ds/4.0)*(K1K+3.0*K2K);
  		Dee.noalias() += (ds/4.0)*(K1D+3.0*K2D);
	}

	Kee = ((Kee.array().abs() > TOL*Kee.norm()).select(Kee.array(),0.0));
	Dee = ((Dee.array().abs() > TOL*Dee.norm()).select(Dee.array(),0.0));

	debug("Assembly of global matrices");
}

//---------------------------------------------------
//------------------------------ build gloabl system
//---------------------------------------------------
void Cosserat::buildLagrangianModel(
	Vnd x1, 
	Vnd x2)
{
	Mad K1J(6,Na*Nm), K2J(6,Na*Nm);
	Mad K1Jt(6,Na*Nm), K2Jt(6,Na*Nm);
	Mnd K1M(Na*Nm,Na*Nm), K2M(Na*Nm,Na*Nm);
	Mnd K1Mt(Na*Nm,Na*Nm), K2Mt(Na*Nm,Na*Nm);
	Mnd K1C(Na*Nm,Na*Nm), K2C(Na*Nm,Na*Nm);
	Mad_ K1Me(Na*Nm,6), K2Me(Na*Nm,6);
	Mad_ K1Ce(Na*Nm,6), K2Ce(Na*Nm,6);
	Vnd K1G(Na*Nm), K2G(Na*Nm);
	V6d tmp;
	M6d AdgInv;

	double K1Vg, K2Vg;

	V13d K1, K2;
	V13d x, dx;

	// initial matrices
	Jb.setZero();
	Jbt.setZero();
	M.setZero();
	Mt.setZero();
	C.setZero();
	G.setZero();
	Ce.setZero();
	Me.setZero();
	x.setZero();
	Vg = 0.0;
	x(0) = 1.0;

	// set states
	q.noalias()  = x1;
	dq.noalias() = x2;

	// compute forward integration step
	double ds = (1.0*(L))/(1.0*(IntegrationSteps));
	double s  = 0.0;

	// do spatial integration
	for (int i = 0; i < IntegrationSteps; i++){

  		lagrangianODE(s,x,Jb,Jbt,K1,K1J,K1Jt,K1M,K1C,K1G,K1Mt,K1Me,K1Ce,K1Vg);

 		lagrangianODE(s+(2.0/3.0)*ds,
 				x+(2.0/3.0)*ds*K1, 
 				Jb+(2.0/3.0)*ds*K1J,
 				Jbt+(2.0/3.0)*ds*K1Jt,
 				K2,K2J,K2Jt,K2M,K2C,K2G,K2Mt,K2Me,K2Ce,K2Vg);

 		s += ds;
  		x.noalias()   += (ds/4.0)*(K1+3.0*K2); 
  		Jb.noalias()  += (ds/4.0)*(K1J+3.0*K2J); 
  		Jbt.noalias() += (ds/4.0)*(K1Jt+3.0*K2Jt); 
  		M.noalias()   += (ds/4.0)*(K1M+3.0*K2M); 
  		C.noalias()   += (ds/4.0)*(K1C+3.0*K2C); 
  		G.noalias()   += (ds/4.0)*(K1G+3.0*K2G); 
  		Mt.noalias()  += (ds/4.0)*(K1Mt+3.0*K2Mt); 
  		Me.noalias()  += (ds/4.0)*(K1Me+3.0*K2Me); 
  		Ce.noalias()  += (ds/4.0)*(K1Ce+3.0*K2Ce); 
  		Vg += (ds/4.0)*(K1Vg+3.0*K2Vg);
	}

	// return configuration and velocities
	g.noalias()   = x.block(0,0,7,1);
	eta.noalias() = x.block(7,0,6,1);

	// compute adjoint actions
	AdmapInv(g,AdgInv);

	// transform Jacobian to local frame
	K1J.noalias() = AdgInv*Jb;
	Jb.noalias() = K1J;

	K1Jt.noalias() = AdgInv*Jbt;
	Jbt.noalias() = K1Jt;
}

//---------------------------------------------------
//----------------------------- build global system
//---------------------------------------------------
void Cosserat::buildInertiaMatrix(
	Vnd v
)
{
	Mad K1J(6,Na*Nm), K2J(6,Na*Nm);
	Mnd K1M(Na*Nm,Na*Nm), K2M(Na*Nm,Na*Nm);

	V7d K1, K2;
	V7d x, dx;

	// initial matrices
	Jb.setZero();
	M.setZero();
	x.setZero();
	x(0) = 1.0;

	// set states
	q.noalias()  = v;

	// compute forward integration step
	double ds = (1.0*(L))/(1.0*(IntegrationSteps));
	double s  = 0.0;

	// do spatial integration
	for (int i = 0; i < IntegrationSteps; i++){

		inertiaODE(s,x,Jb,K1, K1J,K1M);

 		inertiaODE(s+(2.0/3.0)*ds,
 				x+(2.0/3.0)*ds*K1, 
 				Jb+(2.0/3.0)*ds*K1J,
 				K2,K2J,K2M);

 		s += ds;
  		x.noalias()  += (ds/4.0)*(K1+3.0*K2);  
  		Jb.noalias() += (ds/4.0)*(K1J+3.0*K2J);
  		M.noalias()  += (ds/4.0)*(K1M+3.0*K2M);
	}
}

//---------------------------------------------------
//-------------------------- build global system ODE
//---------------------------------------------------
void Cosserat::globalSystemMatODE(
	double s, 
	Mnd &K, 
	Mnd &D)
{
	MatrixXd PMat(Na, Na*Nm);

	// evaluate strain-field
	Phi.eval(s,PMat);

	K.noalias() = ((Ba*PMat).transpose())*Ktt*(Ba*PMat);
	D.noalias() = ((Ba*PMat).transpose())*Dtt*(Ba*PMat);
}

//---------------------------------------------------
//---------------------- build lagrangian system ODE
//---------------------------------------------------
void Cosserat::lagrangianODE(
	double s, 
	V13d x, 
	Mad J, 
	Mad Jt,
	V13d &dx, 
	Mad &dJ, 
	Mad &dJt, 
	Mnd &dM, 
	Mnd &dC, 
	Vnd &dG,
	Mnd &dMt, 
	Mad_ &dMe, 
	Mad_ &dCe,
	double &dVg)
{
	// initialize Shape-function matrix
	MatrixXd PMat(Na, Na*Nm);

	// temporary variable matrix
	Mnd X0(6,Na*Nm);
	M6d X1, X2;
	M4d A;
	M3d R;
	V4d quat;
	V6d xi, eta, r;

	//evaluate strain-field
	Phi.eval(s,PMat);

	// decomposition of configuration space
	xi.noalias()   = (Ba*PMat)*q + Xi0;
	g.noalias()    = x.block(0,0,7,1);
	quat.noalias() = x.block(0,0,4,1);
	eta.noalias()  = x.block(7,0,6,1);

	SE3pos(g, r);
	quat2rot(quat,R);
	strainMapping(R*xi.block(0,0,3,1),A);

	// precompute temporaries
	admap(xi,X2);

	// compute local D-configuration
	(dx.block(0,0,4,1)).noalias() = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias() = R*xi.block(3,0,3,1);
	(dx.block(7,0,6,1)).noalias() = -X2*eta + (Ba*PMat)*dq;

	// overwrite temporary matrices
	Admap(g,X1);
	admap(eta,X2);

	//compute local D-Jacobian
	dJ.noalias()  = X1*(Ba*PMat);
	dJt.noalias() = X1*X2*(Ba*PMat);

	// overwrite Ad(g) -> Ad(g^-1)
	AdmapInv(g,X1);

	X0.noalias() = X1*J;

	// compute local D-inertia matrix
	dM.noalias() = (X0).transpose()*Mtt*(X0);

	// compute local C-coriolis matrix
	//dC.noalias() = (X0).transpose()*((Mtt*X2 - 
	// 	X2.transpose()*Mtt)*(X0) + Mtt*(X0));
	dC.noalias() = (X0).transpose()*((Mtt*X2 - 
	 	X2.transpose()*Mtt)*(X0) + Mtt*(X1*Jt));

	// compute local G-potential vector
	dG.noalias() = (X0).transpose()*Mtt*X1*gvec;

	// compute local time-derivative of M-inertia matrix
	dMt.noalias() = (X1*Jt).transpose()*Mtt*(X0) + (X0).transpose()*Mtt*(X1*Jt);

	// compute local moving-base at \eta_0: Me and Ce matrix
	dMe.noalias() = (X0).transpose()*Mtt;
	dCe.noalias() = (X0).transpose()*(Mtt*X2 - X2.transpose()*Mtt);

	dVg = A11*Rho*(r).dot(gvec);
}

//---------------------------------------------------
//-------------------------- build global system ODE
//---------------------------------------------------
void Cosserat::inertiaODE(
	double s, 
	V7d x, 
	Mad J, 
	V7d &dx, 
	Mad &dJ, 
	Mnd &dM)
{
	// initialize Shape-function matrix
	MatrixXd PMat(Na, Na*Nm);

	// temporary variable matrix
	M6d X1, X2;
	M4d A;
	M3d R;
	V4d quat;
	V6d xi;

	//evaluate strain-field
	Phi.eval(s,PMat);

	// decomposition of configuration space
	xi.noalias()   = (Ba*PMat)*q + Xi0;
	g.noalias()    = x.block(0,0,7,1);
	quat.noalias() = x.block(0,0,4,1);

	quat2rot(quat,R);
	strainMapping(R*xi.block(0,0,3,1),A);

	// precompute temporaries
	admap(xi,X2);

	// compute local D-configuration
	(dx.block(0,0,4,1)).noalias() = (1.0/(2*quat.norm()))*A*quat;
	(dx.block(4,0,3,1)).noalias() = R*xi.block(3,0,3,1);

	// overwrite temporary matrices
	Admap(g,X1);

	//compute local D-Jacobian
	dJ.noalias() = X1*(Ba*PMat);

	// overwrite Ad(g) -> Ad(g^-1)
	AdmapInv(g,X1);

	// compute local D-inertia matrix
	dM.noalias() = (X1*J).transpose()*Mtt*(X1*J);
}
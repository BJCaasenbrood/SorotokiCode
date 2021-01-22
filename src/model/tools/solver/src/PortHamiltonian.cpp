#include "PortHamiltonian.h"

//--------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
PortHamiltonian::PortHamiltonian(
	const char* str)
{

	ConfigFile cf(str);	

  	Na = static_cast<int>(cf.Value("model","NDOF"));
  	Nm = static_cast<int>(cf.Value("model","NMODE"));

  	EC1 = static_cast<bool>(cf.Value("options","ENERGY_CONTROLLER"));

  	h = cf.Value("solver","TIMESTEP");
  	T = cf.Value("solver","TDOMAIN");

	q.noalias() = VectorXd::Zero(Na*Nm);
  	p.noalias() = VectorXd::Zero(Na*Nm);
  	dHdq.noalias() = VectorXd::Zero(Na*Nm);
  	dHdp.noalias() = VectorXd::Zero(Na*Nm);

	q.noalias() = readVecXf("log/state.txt");
	p.noalias() = readVecXf("log/momenta.txt");

	debug("setup PortHamiltonian class");

  	rod.setup(str);
  	control.setup(str);

  	control.Sa.noalias() = rod.Sa;
	
	t = 0.0;
}

//---------------------------------------------------
//----------------- evaluate gradient of Hamiltonian
//---------------------------------------------------
void PortHamiltonian::SolverUpdate()
{
	// assign new port-variables
	Vnd q_(Na*Nm), p_(Na*Nm);

	#if true
		ImplicitTrapzoidStep(q,p,q_,p_);
	#else
		//StromerVerletUpdateStep(q,p,q_,p_);
	#endif

	t += h;
	q.noalias() = q_;
	p.noalias() = p_;
	dHdq.noalias() = rod.gradHq;
  	dHdp.noalias() = rod.gradHp;
}

//---------------------------------------------------
//----------------- evaluate gradient of Hamiltonian
//---------------------------------------------------
void PortHamiltonian::StromerVerletUpdateStep(
	Vnd Q0,
	Vnd P0,
	Vnd &Q1,
	Vnd &P1)
{
	Vnd P_(Na*Nm);

	// evaluate n-point PH-model 
	rod.build(Q0,P0);

	//mid-point backward momenta-step
	P_.noalias() = P0 - 0.5*h*rod.gradHq;

	// full-point forward state-step
	Q1.noalias() = Q0 + h*rod.M.householderQr().solve(P_);

	// evaluate n+1-point PH-model 
	rod.build(Q1,P_);
	P1.noalias() = P_ - 0.5*h*rod.gradHq;
}

//---------------------------------------------------
//----------------- evaluate gradient of Hamiltonian
//---------------------------------------------------
void PortHamiltonian::ImplicitTrapzoidStep(
	Vnd Q0,
	Vnd P0,
	Vnd &Q1,
	Vnd &P1)
{
	Vnd R(2*Na*Nm), dr(2*Na*Nm);

	// evaluate n-point PH/LAG-model
	rod.build(Q0,P0);

	control.q.noalias() 	= Q0;
	control.p.noalias() 	= P0;
	control.t				= t;
	control.dHdq.noalias()  = rod.gradHq;
	control.dHdp.noalias()  = rod.gradHp;
	control.Dee.noalias()	= rod.Dee;

	control.EnergyShaping(rod.Jb,rod.g);

	// computer residual
	R.block(0,0,Na*Nm,1).noalias()     = +h*rod.gradHp;
	R.block(Na*Nm,0,Na*Nm,1).noalias() = -h*rod.gradHq;
	R.block(Na*Nm,0,Na*Nm,1).noalias() += -h*rod.Dee*rod.gradHp;

	if (EC1){
		R.block(Na*Nm,0,Na*Nm,1).noalias() += +h*rod.Sa.transpose()*control.u;
	}

	// hessian update 
	HessianInverse(-R,dr);

	Q1.noalias() = Q0 - dr.block(0,0,Na*Nm,1);
	P1.noalias() = P0 - dr.block(Na*Nm,0,Na*Nm,1);
}

//---------------------------------------------------
//--------------------------- compute Hessian matrix
//---------------------------------------------------
void PortHamiltonian::HessianInverse(
	Vnd R, 
	Vnd &dr)
{

	int n = Na*Nm;
	double alpha = (2.0/3.0);
	double beta = 1.0;

	if (EC1){
		beta = control.Kd;
	}

	Mnd S(2*n,2*n);
	Mnd Mi(n,n);

	Mi.noalias() = rod.M.partialPivLu().solve(Mnd::Identity(n,n));

	S.block(0,n,n,n).noalias() = -alpha*h*Mi;
	S.block(0,0,n,n).noalias() = Mnd::Zero(n,n);
	S.block(n,0,n,n).noalias() = alpha*h*rod.Kee;
	S.block(n,n,n,n).noalias() = alpha*h*(1+beta)*rod.Dee*Mi;
	S.block(n,n,n,n).noalias() += -alpha*h*(rod.Mt - rod.C)*Mi;
	S.noalias() += Mnd::Identity(2*n,2*n);

	dr.noalias() = S.partialPivLu().solve(R);
}
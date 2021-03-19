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
  	Ny = 1;
  	Nu = 1;

  	EC1 = static_cast<bool>(cf.Value("options","ENERGY_CONTROLLER"));

  	h = cf.Value("solver","TIMESTEP");
  	T = cf.Value("solver","TDOMAIN");

	q.noalias() = VectorXd::Zero(Na*Nm);
  	p.noalias() = VectorXd::Zero(Na*Nm);
  	q_.noalias() = VectorXd::Zero(Na*Nm);
  	p_.noalias() = VectorXd::Zero(Na*Nm);
  	z.noalias() = VectorXd::Zero(7);
  	z_.noalias() = VectorXd::Zero(6);
  	//q_(0) = 1e-2;

  	dHdq.noalias() = VectorXd::Zero(Na*Nm);
  	dHdp.noalias() = VectorXd::Zero(Na*Nm);
  	A.noalias() = MatrixXd::Zero(2*Na*Nm,2*Na*Nm);
  	B.noalias() = MatrixXd::Zero(2*Na*Nm,Nu);
  	C.noalias() = MatrixXd::Zero(Ny,2*Na*Nm);

	q.noalias() = readVecXf("log/state.txt");
	p.noalias() = readVecXf("log/momenta.txt");
	q_.noalias() = readVecXf("log/estimate.txt");

	debug("setup PortHamiltonian class");

  	rod.setup(str);
  	control.setup(str);

  	control.Sa.noalias() = rod.Sa;
}

//---------------------------------------------------
//----------------- evaluate gradient of Hamiltonian
//---------------------------------------------------
void PortHamiltonian::SolverUpdate()
{
	// assign new port-variables
	Vnd qn(Na*Nm), pn(Na*Nm);

	#if true
		ImplicitTrapzoidStep(q,p,qn,pn);
	#else
		StromerVerletUpdateStep(q,p,qn,pn);
	#endif

	q.noalias() = qn;
	p.noalias() = pn;
	dHdq.noalias() = rod.gradHq;
  	dHdp.noalias() = rod.gradHp;
  	z.noalias()  = rod.g;
  	z_.noalias() = rod.eta;
}

//---------------------------------------------------
//----------------- evaluate gradient of Hamiltonian
//---------------------------------------------------
void PortHamiltonian::TimeUpdate(){
	t +=h;
}

//---------------------------------------------------
//----------------- evaluate gradient of Hamiltonian
//---------------------------------------------------
void PortHamiltonian::ObserverUpdate()
{
	// assign new port-variables
	Vnd qn(Na*Nm), pn(Na*Nm);

	#if true
		LinearizedStep(q_,p_,qn,pn);
	#else
		//StromerVerletUpdateStep(q,p,qn,pn);
	#endif

	q_.noalias() = qn;
	p_.noalias() = pn;
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
	int n = Na*Nm;
	Vnd Rs(2*n), dr(2*n);

	// evaluate n-point PH/LAG-model
	rod.build(Q0,P0);

	control.t				= t;
	control.q.noalias() 	= Q0;
	control.p.noalias() 	= P0;
	control.dHdq.noalias()  = rod.gradHq;
	control.dHdp.noalias()  = rod.gradHp;
	control.Dee.noalias()	= rod.Dee;

	control.EnergyShaping(rod.Jb,rod.g);

	// computer residual
	Rs.block(0,0,n,1).noalias() = +h*rod.gradHp;
	Rs.block(n,0,n,1).noalias() = -h*rod.gradHq;
	Rs.block(n,0,n,1).noalias() += -h*rod.Dee*rod.gradHp;

	if (EC1){
		double alpha = 1;//smoothstep(t,0,1); // prevents fast accelerations.
		Rs.block(n,0,n,1).noalias() += +alpha*h*rod.Sa.transpose()*control.u;
	}

	// hessian update 
	HessianInverse(-Rs,dr);

	Q1.noalias() = Q0 - dr.block(0,0,n,1);
	P1.noalias() = P0 - dr.block(n,0,n,1);

	// compute joint accelerations
	rod.ddq.noalias() = rod.M.partialPivLu().solve(-(1.0f/h)*dr.block(Na*Nm,0,Na*Nm,1) + rod.Mt*rod.dq);

	//Mxd P(2*Na*Nm,2*Na*Nm);
	//Mxd Q(2*Na*Nm,2*Na*Nm);
	//Mxd R(1,1);

	//R << 1.0;
	//Q.noalias() = Mnd::Identity(2*Na*Nm,2*Na*Nm);

	//solveRiccatiArimotoPotter(A,C,Q,R,P);

	//cout << A << endl;
	//cout << C << endl;
	//cout << P << endl;

	//cout << "check: " << A.transpose()*P + P*A  - P*C*C.transpose()*P + Q << endl;

	//sleep(1200);

	//SchurComplement(A.transpose(),C.transpose());
}

//---------------------------------------------------
//----------------- evaluate gradient of Hamiltonian
//---------------------------------------------------
void PortHamiltonian::LinearizedStep(
	Vnd Q0,
	Vnd P0,
	Vnd &Q1,
	Vnd &P1)
{
	
	int n = Na*Nm;
	Vnd Rs(2*n), dr(2*n);

	// evaluate n-point
	rod.build(Q0,P0);

	Mnd Mi(n,n);
	Mi.noalias() = rod.M.partialPivLu().solve(Mnd::Identity(n,n));

	// computer residual
	Rs.block(0,0,n,1).noalias() = +h*rod.gradHp;
	Rs.block(n,0,n,1).noalias() = -h*rod.gradHq;
	Rs.block(n,0,n,1).noalias() += -h*rod.Dee*rod.gradHp;

	if (EC1){
		double alpha = 1.0;
		Rs.block(n,0,n,1).noalias() += +alpha*h*rod.Sa.transpose()*control.u;
	}

	control.KalmanCorrection(rod.Jb,rod.g,z,rod.eta - z_);

	Rs.block(n,0,n,1).noalias() += +h*control.uk;

	// hessian update 
	HessianInverse(-Rs,dr);

	Q1.noalias() = Q0 - dr.block(0,0,n,1);
	P1.noalias() = P0 - dr.block(n,0,n,1);

}


//---------------------------------------------------
//--------------------------- compute Hessian matrix
//---------------------------------------------------
void PortHamiltonian::HessianInverse(
	Vnd Rs, 
	Vnd &dr){

	int n = Na*Nm;
	double alpha = (2.0/3.0);
	double beta = 1.0;

	if (EC1){
		beta = control.Kd;
	}

	Mnd Mi(n,n);
	Mi.noalias() = rod.M.partialPivLu().solve(Mnd::Identity(n,n));

	A.block(0,n,n,n).noalias() = Mi;
	A.block(0,0,n,n).noalias() = Mnd::Zero(n,n);
	A.block(n,0,n,n).noalias() = -rod.Kee;
	A.block(n,n,n,n).noalias() = -(1+beta)*rod.Dee*Mi;
	A.block(n,n,n,n).noalias() += (rod.Mt - rod.C)*Mi;

	dr.noalias() = (-alpha*h*A + Mnd::Identity(2*n,2*n)).partialPivLu().solve(Rs);
}
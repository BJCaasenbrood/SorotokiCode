#include "PortController.h"

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
PortController::PortController()
{

}

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
void PortController::setup(
	const char* str)
{

	ConfigFile cf(str);	

	double Xd, Yd, Zd;
  	double Q1d, Q2d, Q3d, Q4d;
  	int Order;

  	Na = static_cast<int>(cf.Value("model","NDOF"));
  	Nm = static_cast<int>(cf.Value("model","NMODE"));
  	h = cf.Value("solver","TIMESTEP");
  	T = cf.Value("solver","TDOMAIN");

  	Lambda = cf.Value("control","LAMBDA");
  	LambdaK = cf.Value("control","LAMBDAK");
  	Order = cf.Value("control","SPLINEORDER");
  	Kp = cf.Value("control","KP");
  	Kd = cf.Value("control","KD");
  	L  = cf.Value("control","LK");
  	K1 = cf.Value("control","KF1");
  	K2 = cf.Value("control","KF2");
  	Xd = cf.Value("setpoint","Xd");
  	Yd = cf.Value("setpoint","Yd");
  	Zd = cf.Value("setpoint","Zd");
  	Q1d = cf.Value("setpoint","Q1d");
  	Q2d = cf.Value("setpoint","Q2d");
  	Q3d = cf.Value("setpoint","Q3d");
  	Q4d = cf.Value("setpoint","Q4d");

	q.noalias()  = VectorXd::Zero(Na*Nm);
  	p.noalias()  = VectorXd::Zero(Na*Nm);
  	u.noalias()  = VectorXd::Zero(Na*Nm);
  	uk.noalias() = VectorXd::Zero(Na*Nm);
  	dHdq.noalias() = VectorXd::Zero(Na*Nm);
  	dHdp.noalias() = VectorXd::Zero(Na*Nm);
  	H.noalias()  = VectorXd::Zero(2);

  	Mnd MatPolyX = readMatXf("log/splineXd.txt");

  	VectorXd Xs(MatPolyX.rows());
  	VectorXd Ys(MatPolyX.rows());

  	Xs.noalias() = MatPolyX.block(0,0,MatPolyX.rows(),1);
  	Ys.noalias() = MatPolyX.block(0,1,MatPolyX.rows(),1);

  	const auto fit = SplineFitting<S1d>::Interpolate(Ys.transpose(), 
	  	1,//std::min<int>(Order, 3),
    	scaled_values(Xs).transpose());

  	S1d S(fit);

 	spline_X = S;

  	gd << Q1d, Q2d, Q3d, Q4d, Xd, Yd, Zd;

	debug("setup PortController class");
	
	t = 0.0;
}

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
void PortController::EnergyShaping(
	Mad J,
	V7d g){

	int n = Na*Nm;
	int m = n;

	M6d K;
	V6d R;
	V6d k;

	Mnd N(n,n);
	Vnd dx(n);
	M4d G, Gi;
	M6d Tang;
	V6d Phi;

	R.setZero();
	K.setZero();
	dx.setZero();

	// set virtual stiffness gains
	k << K1,K1,K1,K2,K2,K2;
	K.diagonal() = k;

	// relative spatial-twist between g(l) and g*
	SE3(g,G);
	SE3Inv(gd,Gi);
	logmapSE3(Gi*G,Phi);
	//tmapSE3(Phi,Tang);

	// compute residual spring force
	//R.noalias() = K*Tang*Phi*smoothstep(t,0.0,2);
	R.noalias() = K*Phi*smoothstep(t,0.0,2);

	//compute dq - desired potential energy 
	dx.noalias() = J.transpose()*(J*J.transpose() + 
		Lambda*Mnd::Identity(6,6)).householderQr().solve(R);

	// energy-shaping + damping injection -- - Kp*dx 
	u.noalias() = Sa*(dHdq - Kp*dx - Kd*Dee*dHdp);

	H(0) = 0;
	H(0) += dx.transpose()*Kp*q;
	//H(0) += dHdq.transpose()*q;
	//H(0) += dHdq.transpose()*q;
}                         

//---------------------------------------------------
//-------------------------------- class constructor
//---------------------------------------------------
void PortController::KalmanCorrection(
	Mad J,
	V7d g,
	V7d z,
	V6d ed){

	int n = Na*Nm;
	int m = n;

	M6d K;
	V6d R;
	V6d k;

	Mnd N(n,n);
	Vnd dx(n);
	M4d G, Gi;
	M6d Tang;
	V6d Phi;

	R.setZero();
	K.setZero();
	dx.setZero();

	// set virtual stiffness gains
	k << K1,K1,K1,K2,K2,K2;
	K.diagonal() = k;

	// relative spatial-twist between g(l) and measurement
	SE3(g,G);
	SE3Inv(z,Gi);
	logmapSE3(Gi*G,Phi);
	tmapSE3(Phi,Tang);

	// compute residual spring force
	R.noalias() = K*Tang*Phi;

	//compute dq - desired potential energy 
	//x.noalias() = J.transpose()*(J*J.transpose() + 
	//	LambdaK*Mnd::Identity(6,6)).householderQr().solve(R);
	dx.noalias() = J.transpose()*R;

	// kalman output
	uk.noalias() = -L*dx - LambdaK*J.transpose()*ed;
}                         
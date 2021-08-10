[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 1 
WRITE_OUTPUT         = 1 
ACTUSPACE =0
PROPERTYSET =-1
DISCONTINIOUS        = 1 

[cosserat] 
K1 = 0
K2 = 1
K3 = 1
E1 = 0
E2 = 0
E3 = 0

[model] 
NMODE   = 6
NDISC   = 2
NDOF    = 2

[solver] 
TDOMAIN = 25
SPACESTEP = 50
TIMESTEP  = 0.0066667
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
MAX_IMPL  = -1
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80
ADAMPING  = 1

[physics] 
LENGTH   = 0.36
RHO      = 550
EMOD     = 50
NU       = 0.4
MU       = 0.45
PRS_AREA = 1e-5
GRAVITY  = 0  0  0
RADIUS   = 0.01
AREA     = 0.007854
J_XX     = 0.00125
J_YY     = 0.000625
J_ZZ     = 0.000625

[control] 
KP = 100
KD = 1
LK = 0
KF1 = 0
KF2 = 1
LAMBDA    = 10
LAMBDAK   = 0
SPLINEORDER    = 1

[setpoint] 
Q1d = 0
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  0.2
Yd =  0.10806
Zd =  0.16829

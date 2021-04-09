[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 1 
WRITE_OUTPUT         = 1 
ACTUSPACE =1
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
NMODE   = 8
NDISC   = 2
NDOF    = 2

[solver] 
TDOMAIN = 15
SPACESTEP = 150
TIMESTEP  = 0.02
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
MAX_IMPL  = -1
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80
ADAMPING  = 1

[physics] 
LENGTH   = 0.12
RHO      = 1200
EMOD     = 25
NU       = 0.4
MU       = 0.1
PRS_AREA = 1e-5
GRAVITY  = 0           0       -9.81
RADIUS   = 0.01
AREA     = 0.00031416
J_XX     = 5e-05
J_YY     = 2.5e-05
J_ZZ     = 2.5e-05

[control] 
KP = 0.1
KD = 15
LK = 0
KF1 = 1e-12
KF2 = 1
LAMBDA    = 0.002
LAMBDAK   = 0
SPLINEORDER    = 1

[setpoint] 
Q1d = 1
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  0.05
Yd =  0
Zd =  0.05
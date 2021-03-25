[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 1 
WRITE_OUTPUT         = 1 
ACTUSPACE =1
PROPERTYSET =-1
DISCONTINIOUS        = 0 

[cosserat] 
K1 = 0
K2 = 1
K3 = 0
E1 = 1
E2 = 0
E3 = 0

[model] 
NMODE   = 4
NDISC   = 1
NDOF    = 2

[solver] 
TDOMAIN = 15
SPACESTEP = 10
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
LENGTH   = 0.063
RHO      = 50
EMOD     = 250
NU       = 0.4
MU       = 0.05
PRS_AREA = 1e-5
GRAVITY  = 0       -9.81           0
RADIUS   = 0.01
AREA     = 0.0028274
J_XX     = 0.00045
J_YY     = 0.000225
J_ZZ     = 0.000225

[control] 
KP = 0.5
KD = 0
LK = 0.5
KF1 = 1e-09
KF2 = 1
LAMBDA    = 0.005
LAMBDAK   = 0.005
SPLINEORDER    = 1

[setpoint] 
Q1d = 0
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  0.060411
Yd =  0
Zd =  -0.038467

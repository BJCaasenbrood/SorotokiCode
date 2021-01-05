[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 1 
WRITE_OUTPUT         = 1 
ACTUATION_SPACE =1
DISCONTINIOUS        = 0 

[cosserat] 
K1 = 0
K2 = 1
K3 = 0
E1 = 1
E2 = 0
E3 = 0

[model] 
NMODE   = 2
NDISC   = 1
SDOMAIN = 0.063
TDOMAIN = 15

[solver] 
SPACESTEP = 5
TIMESTEP  = 0.033333
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
LAMBDA    = 1
MAX_IMPL  = -1
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80
ADAMPING  = 1

[physics] 
RHO      = 10
EMOD     = 50
NU       = 0.4
MU       = 0.3
PRS_AREA = 1e-5
GRAVITY  = -9.81           0           0
RADIUS   = 0.01
AREA     = 0.00031416
J_XX     = 5e-07
J_YY     = 2.5e-07
J_ZZ     = 2.5e-07

[control] 
KP = 1
KD = 0.1

[setpoint] 
Q1d = 1
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  0.058101
Yd =  0
Zd =  -0.03398

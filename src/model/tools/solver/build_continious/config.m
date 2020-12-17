[options] 
KINEMATIC_CONTROLLER = 1 
ENERGY_CONTROLLER    = 0 
WRITE_OUTPUT         = 1 
ACTUATION_SPACE =-1
DISCONTINIOUS        = 0 

[cosserat] 
K1 = 1
K2 = 1
K3 = 0
E1 = 0
E2 = 0
E3 = 0

[model] 
NMODE   = 5
NDISC   = 1
SDOMAIN = 1
TDOMAIN = 40

[solver] 
SPACESTEP = 12
TIMESTEP  = 0.001
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
LAMBDA    = 0.001
MAX_IMPL  = -1
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80
ADAMPING  = 1

[physics] 
RHO      = 0.01
EMOD     = 50
NU       = 0.4
MU       = 4
PRS_AREA = 1e-5
GRAVITY  = 0
RADIUS   = 0.01
AREA     = 0.00031416
J_XX     = 1.3333e-06
J_YY     = 3.3333e-07
J_ZZ     = 0.00014167

[control] 
KP = 0.01
KD = 0.05

[setpoint] 
Q1d = 0
Q2d = 0
Q3d = 0.95106
Q4d = -0.30902
Xd =  0.5
Yd =  0
Zd =  -0.1

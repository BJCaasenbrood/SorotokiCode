[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 0 
WRITE_OUTPUT         = 1 
ACTUATION_SPACE =-1
DISCONTINIOUS        = 0 

[cosserat] 
K1 = 0
K2 = 1
K3 = 0
E1 = 0
E2 = 0
E3 = 0

[model] 
NMODE   = 8
NDISC   = 1
NDOF   = 1
SDOMAIN = 0.12
TDOMAIN = 20

[solver] 
TDOMAIN = 20
SPACESTEP = 25
TIMESTEP  = 0.1
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
LAMBDA    = 2
MAX_IMPL  = 
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80
ADAMPING  = 1

[physics] 
LENGTH   = 0.12
RHO      = 150
EMOD     = 200
NU       = 0.4
MU       = 0.2
PRS_AREA = 1e-5
GRAVITY  = 0           0        9.81
RADIUS   = 0.01
AREA     = 0.00031416
J_XX     = 5e-07
J_YY     = 2.5e-07
J_ZZ     = 2.5e-07

[control] 
KP = 1
KD = 0

[setpoint] 
Q1d = 1
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  0.05
Yd =  0
Zd =  -0.01

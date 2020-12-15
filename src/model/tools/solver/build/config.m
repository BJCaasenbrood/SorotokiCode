[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 1 
WRITE_OUTPUT         = 1 
ACTUATION_SPACE =1
DISCONTINIOUS        = 1 

[cosserat] 
K1 = 0
K2 = 0
K3 = 1
E1 = 0
E2 = 0
E3 = 0

[model] 
NMODE   = 5
NDISC   = 2
NDOF   = 1
SDOMAIN = 0.36
TDOMAIN = 15

[solver] 
TDOMAIN = 15
SPACESTEP = 30
TIMESTEP  = 0.083333
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
LAMBDA    = 0.02
MAX_IMPL  = 
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80
ADAMPING  = 1

[physics] 
LENGTH   = 0.36
RHO      = 10
EMOD     = 250
NU       = 0.4
MU       = 0.83
PRS_AREA = 1e-5
GRAVITY  = -9.81           0           0
RADIUS   = 0.01
AREA     = 0.00031416
J_XX     = 5e-07
J_YY     = 2.5e-07
J_ZZ     = 2.5e-07

[control] 
KP = 0.001
KD = 0.1

[setpoint] 
Q1d = 1
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  0.15
Yd =  0.05
Zd =  -0.15

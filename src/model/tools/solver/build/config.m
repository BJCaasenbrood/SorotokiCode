[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 0 
WRITE_OUTPUT         = 1 
DISCONTINIOUS        = 0 

[cosserat] 
K1 = 0
K2 = 1
K3 = 1
E1 = 0
E2 = 0
E3 = 0

[model] 
NMODE   = 1
NDISC   = 1
SDOMAIN = 1
TDOMAIN = 2

[solver] 
SPACESTEP = 11
TIMESTEP  = 0.0001
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
LAMBDA    = 1e-2
MAX_IMPL  = -1
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80

[physics] 
RHO      = 0.7656
EMOD     = 2000
NU       = 0.4
MU       = 0.005
PRS_AREA = 1e-5
GRAVITY  = 0
RADIUS   = 0.01
AREA     = 0.0012566
J_XX     = 0.00015312
J_YY     = 7.656e-05
J_ZZ     = 7.656e-05

[control] 
KP = 1e-07
KD = 0.15

[setpoint] 
Q1d = 1
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  1
Yd =  0
Zd =  0

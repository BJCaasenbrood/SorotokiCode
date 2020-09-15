[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 1 
WRITE_OUTPUT         = 1 
DISCONTINIOUS        = 1 

[cosserat] 
K1 = 0
K2 = 0
K3 = 1
E1 = 0
E2 = 0
E3 = 0

[model] 
NMODE   = 4
NDISC   = 2
SDOMAIN = 1
TDOMAIN = 25

[solver] 
SPACESTEP = 21
TIMESTEP  = 0.03
INTSTEP   = 100
ATOL      = 1e-2
RTOL      = 1e-2
LAMBDA    = 1e-2
MAX_IMPL  = -1
MAX_ITER  = 1
MAX_IK    = 1500
SPEEDUP   = 80

[physics] 
RHO      = 0.01
EMOD     = 1000000
NU       = 0.4
MU       = 0.5
PRS_AREA = 1e-5
GRAVITY  = -9.81
RADIUS   = 0.01
AREA     = 3.0000e-04
J_XX     = 2.7250e-10
J_YY     = 2.2500e-11
J_ZZ     = 2.5000e-10

[control] 
KP = 0.0005
KD = 0

[setpoint] 
Q1d = 1
Q2d = 0
Q3d = 0
Q4d = 0
Xd =  0.8
Yd =  0
Zd =  0

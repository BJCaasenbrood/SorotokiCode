[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 1 
WRITE_OUTPUT         = 1 
DISCONTINIOUS        = 1 

[cosserat] 
K1 = 0
K2 = 1
K3 = 0
E1 = 0
E2 = 0
E3 = 0

[model] 
NMODE   = 4
NDISC   = 3
SDOMAIN = 1
TDOMAIN = 25

[solver] 
SPACESTEP = 21
TIMESTEP  = 0.01
INTSTEP   = 100
ATOL      = 1e-2
LAMBDA    = 1e-2
MAX_IMPL  = -1
MAX_IK    = 1500
SPEEDUP   = 80

[physics] 
RHO      = 21
EMOD     = 0.01
NU       = 100
MU       = 1e-2
PRS_AREA = 1e-2
GRAVITY  = -1
RADIUS   = 1500
AREA     = 3.0000e-04
J_XX     = 2.7250e-10
J_YY     = 2.2500e-11
J_ZZ     = 2.5000e-10

[control] 
KP = 2.2500e-11
KD = 2.5000e-10

[setpoint] 
Xd =  2.2500e-11
Yd =  2.2500e-11
Q1d = 2.2500e-11
Q1d = 2.2500e-11
Q1d = 
Q1d = 0.0

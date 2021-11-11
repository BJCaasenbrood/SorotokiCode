[options] 
KINEMATIC_CONTROLLER = 0 
ENERGY_CONTROLLER    = 0 
WRITE_OUTPUT         = 1 
ACTUSPACE =1
PROPERTYSET =-1
DISCONTINIOUS        = 0 

[cosserat] 
K1 = 0
K2 = 1
K3 = 1
E1 = 1
E2 = 0
E3 = 0

[model] 
NMODE   = 1
NDISC   = 1
NDOF    = 3

[solver] 
TDOMAIN = 3
SPACESTEP = 25
TIMESTEP  = 0.001
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
RHO      = 1800
EMOD     = 1500
NU       = 0.4
MU       = 0.05
PRS_AREA = 1e-5
GRAVITY  = 0  0  0
RADIUS   = 0.01
AREA     = 0.00031416
J_XX     = 5e-05
J_YY     = 2.5e-05
J_ZZ     = 2.5e-05

[control] 
KP = 0.0005
KD = 20

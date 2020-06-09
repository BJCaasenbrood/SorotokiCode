% simulink settings

t_span = 5;           % simulation time
dt     = 1.00e-03;    % time-steps for solver

% general physical parameters
m      = 0.0494;        % mass of soft robot
g      = 9.81;        % gravitation
r      = 0.0235;       % radius of cross-section
d0     = 0.0644;       % relaxation length
alpha0 = 0;           % relaxation angle (unused)
I      = eye(3,3);    % identity matrix
c1     = 0.12;           % damping for elongation
c2     = 5.93e-7;       % damping for curvature kx
c3     = 5.93e-7;       % damping for curvature ky

% set noise levels

noise_q     = 0*[1e-6, 1e-1, 1e-1];
noise_dq    = 0*[1e-4, 1e-4, 1e-4];
noise_ddq 	= 0*[1e-4, 1e-4, 1e-4];

% initialize states

dq_0 = [0; 0; 500];   % initial conditions for states
q_0  = [d0; 10; .1];  % initial conditions for state velocities 

% initialize nonlinear stiffness parameters

% ke(e) = b3 - b1*b2*(tanh(b2*e)^2 - 1) : nonlinear elongation stiffness
% kb(b) = b6 - b4*b5*(tanh(b4*e)^2 - 1) : nonlinear bending stiffness

b1 =  2.23e3;         % stiffness parameter: elongation
b2 =  1.74e3;       % stiffness parameter: elongation
b3 =  -4.55e2;       % stiffness parameter: elongation

b4 =  4.23e-1;         % stiffness parameter: bending
b5 =  3.99e-1;         % stiffness parameter: bending
b6 =  -2.29e-1;         % stiffness parameter: bending

% set the parameter vector for stiffness
p_true = [b1; b2; b3; b4; b5; b6];

% exit


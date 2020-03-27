clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',5);

mdl = mdl.set('tspan',10,...
              'Grav',9.81,...
              'Jacobian',false,...
              'MovieAxis',[-1 1 -1 1 -2.5 1]*0.5);

%% generate dynamic model
mdl = mdl.generate();

%% assign controllers
%mdl.Controller = @(t,g) Controller(t,g);
%mdl.Pressure = @(t) [0;0.75;0.75;1.0;0.0;1.0];
mdl.q0 = zeros(5*2,1);
mdl.q0(1) = 1;
mdl.q0(4) = 8;
mdl.q0(7) = 8;

%% simulate soft robot
tic; mdl = mdl.simulate(); toc

%% show simulation
mdl.showModel();

function tau = Controller(t,g)
kp = 1e-4;
tau = zeros(6,1);
tau(4:6) = kp*(g - [0.3;0.5;0.5]*clamp(t/5,1e-2,1));
end
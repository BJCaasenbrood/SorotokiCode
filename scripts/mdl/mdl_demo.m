clr;
%% assign free DOF
tab = [0,1,1,0,0,0];
mdl = Model(tab,'NModal',3);

%% generate dynamic model
mdl = mdl.set('tspan',10,...
              'Grav',0,...
              'MovieAxis',[-1 1 -1 1 -2.5 1]*0.5);

mdl = mdl.generate();

%% assign control law
mdl.Controller = @(t,g) Controller(t,g);

%% simulate soft robot
tic; mdl = mdl.simulate(); toc

%% 
mdl.showModel();

function tau = Controller(t,g)
kp = 1e-4*clamp(t/5,1e-2,1);
kd = 0;


gd = [0.5,0.5,0];%[0.75,0.75*cos(0.5*t),0.75*sin(0.5*t)];

Q = g(1:4);
r = g(5:7);

tau = zeros(6,1);
tau(4:6) = kp*(r' - gd(:));
end
clr;
%% assign free DOF
tab = [0,1,1,0,0,0];
mdl = Model(tab,'NModal',3);

%% generate dynamic model
mdl = mdl.set('tspan',10,'MovieAxis',[-0.5 0.5 -0.5 0.5 -1.25 0.5]);

mdl = mdl.generate();

%% assign control law
mdl.Controller = @(t,g) Controller(t,g);

%% simulate soft robot
tic; mdl = mdl.simulate(); toc

%% 
mdl.showModel();

function F = LoadDistribution(x,t,P1,P2)
t1 = clamp(t/5,0,1);

A = 5e-6;
r = 1;

H = [ 0,              0,             0;
      0, -0.5*r*sqrt(3), 0.5*r*sqrt(3);
     -r,          0.5*r,         0.5*r;
     -0,             -0,            -0;
      0,              0,             0;
      0,              0,             0];

h = 0.05;
a1 = t1;
a2 = t1;
    
F1 = a1*H*P1(:)*(rdelta(x-0.5,h));
F2 = a2*H*P2(:)*(-rdelta(x-0.5,h)+rdelta(x-1+h,h));

F = A*(F1 + F2);
end

function tau = Controller(t,g)
kp = 1e-4*clamp(t/5,0,1);
kd = 0;


gd = [0.75,0.75*cos(0.5*t),0.75*sin(0.5*t)];

Q = g(1:4);
r = g(5:7);

tau = zeros(6,1);
tau(4:6) = kp*(r' - gd(:));
end
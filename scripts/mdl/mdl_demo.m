clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',3);

mdl = mdl.set('tspan',10,...
              'Grav',9.81,...
              'Jacobian',true,...
              'MovieAxis',[-1 1 -1 1 -2.5 1]*0.5);

%% generate dynamic model
mdl = mdl.generate();

%% assign controllers
%mdl.Controller = @(t,g) Controller(t,g);
mdl.Pressure = @(t) 1*[0;.75;0;0;.5;0];

%% simulate soft robot
mdl = mdl.csolve(); 

%% show simulation
figure(12);
plot(mdl.get('t'),mdl.get('g'),'-','linewidth',1.0);

figure(123);
semilogy(mdl.get('t'),vecnorm(gradient(mdl.get('g').')),'-','linewidth',1.0);

figure(121);
plot(mdl.get('t'),mdl.get('tau'),'-','linewidth',1.0);

pause;

mdl.showModel();

function tau = Controller(t,g)
kp = 1e-4;
tau = zeros(6,1);
tau(4:6) = kp*(g - [0.3;0.5;0.5]*clamp(t/5,1e-2,1));
end
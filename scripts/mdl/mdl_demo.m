clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',3);

mdl = mdl.set('tspan',10,...
              'Grav',9.81,...
              'Jacobian',true,...
              'LumpedMass',true,...
              'MovieAxis',[-1 1 -1 1 -2.5 1]*0.5);

%% generate dynamic model
mdl = mdl.generate();

%% assign controllers
%mdl.Controller = @(t,g) Controller(t,g);
mdl.Pressure = @(t) 0.0*[3.0;.0;0;2;2;0];
mdl.q0 = zeros(mdl.NDof*mdl.NModal);
mdl.q0(1) = 5;
mdl.q0(4) = 0;
%% simulate soft robot
mdl = mdl.csolve(); 

%% show simulation
figure(12); hold on;
plot(mdl.get('t'),mdl.get('g'),'-','linewidth',1.0);

figure(123); hold on;
semilogy(mdl.get('t'),vecnorm(abs(gradient(mdl.get('g').'))),'-','linewidth',1.0);
% 
% figure(13); hold on;
% plot(mdl.get('t'),mdl.get('tau'),'-','linewidth',1.0);

pause;
% 
mdl.showModel();

function tau = Controller(t,g)
kp = 1e-4;
tau = zeros(6,1);
tau(4:6) = kp*(g - [0.3;0.5;0.5]*clamp(t/5,1e-2,1));
end
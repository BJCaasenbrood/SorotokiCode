clr;
%% assign free DOF
mdl = Model([0,1,1,1,0,0],'NModal',1);

mdl = mdl.set('tspan',2,...
              'Grav',9.81,...
              'Jacobian',true,...
              'MovieAxis',[-1 1 -1 1 0 1.85]*0.85);

%% generate dynamic model
mdl = mdl.generate();

%% assign controllers
mdl.gain = [3e-4,2e-5];
mdl.point = [];%[0,0,0,.8,0,0.5];
mdl.Pressure = @(t) 0*[2;0;0];
mdl.q0 = zeros(mdl.NDof*mdl.NModal);
mdl.dq0 = zeros(mdl.NDof*mdl.NModal);
mdl.q0(3) = .1;
mdl.dq0(2) = 0;
%% simulate soft robot
mdl = mdl.csolve(); 

%% show simulation
% figure(14); hold on;
% subplot(2,1,1);
t = mdl.get('t');
q = mdl.g;
u = mdl.tau;
ge = mdl.ge;
% plot(t,q,'-','linewidth',1.0);
% 
% subplot(2,1,2);
% plot(t,u,'-','linewidth',1.0);


% figure(15);
% plot3(ge(:,7),ge(:,6),-ge(:,5),'-','linewidth',1.0);
% axis equal;
% t = state(:,1);
% z = state(:,2:end);
% plot(t,z);
% l0 = 0.064;
% g(:,3) = (l0 - q(:,1))/l0;
% g(:,1) = q(:,2)*l0;
% g(:,2) = q(:,3)*l0;
% 
% mdl = mdl.set('t',t);
% mdl = mdl.set('g',g);
% 
mdl.showModel();



function tau = Controller(t,g)
kp = 1e-4;
tau = zeros(6,1);
tau(4:6) = kp*(g - [0.3;0.5;0.5]*clamp(t/5,1e-2,1));
end
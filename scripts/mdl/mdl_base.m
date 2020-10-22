%clr;
%% init. states
l0   = 0.064;               % undeformed length (m)
q0   = [0.064,5,-5];        % initial deformations
dq0  = [0,0,0];             % initial deform. rates
lam0 = [0.0,0.04];               % initial creep strains
p0 = [1.15*2234, 0.95*1740, 1.13*0.42292, 0.85*0.39552];

ShowSimulation = false;
%% ordinary diff. solver (4th order runge kutta)
dt = 1e-1 
tspan = 0:dt:30; 

disp('* simulating dynamics...')
[t,x] = ode15s(@(t,x) DifferentialEq(t,x),tspan,[q0 dq0,lam0,p0]);
%[t,x] = odesolve(@(t,x) DifferentialEq(t,x),tspan,[q0 dq0,lam0,p0]);

%% strip solution matrix x
Elong = (x(:,1) - l0)/l0;
KxDeg = x(:,1).*x(:,2)/(pi/180);
KyDeg = x(:,1).*x(:,3)/(pi/180);
Beta = x(:,1).*x(:,2);%sqrt(x(:,2).^2 + x(:,3).^2).*x(:,1);

y = EndEffector(x(:,1),x(:,2),x(:,3));
%figure;
%set(gcf,'renderer','Painters')
%plot3(y(:,1),y(:,2),y(:,3),'linewidth',2,'color',col(1));

kb = stiffB(Beta,x(:,9:12));
kbt = stiffB(Beta,[2234, 1740, 0.42292, 0.39552]);

ke = stiffE(x(:,1) - l0,x(:,9:12));
ket = stiffE(x(:,1) - l0,[2234, 1740, 0.42292, 0.39552]);

%% show state trajectories

figure(10); subplot(3,1,1);
hold all;
hold on;
plot(t,x(:,1)*1e3,'linewidth',1.5,'color',col(1)); %ylabel('$l - l_0$')
%plot(t,0.07*1000 + 1*sin(t),'--k','linewidth',1.0);
subplot(3,1,2);
hold on;
%plot(t,40*sin(t),'--k','linewidth',1.);
%plot(t,x(:,2),'linewidth',1.5,'color',col(1)); ylabel('$\kappa_x$')
subplot(3,1,3);
hold on; 
%plot(t,40*cos(t),'--k','linewidth',1.);
plot(t,x(:,3),'linewidth',1.5,'color',col(1)); ylabel('$\kappa_y$')
figure(10); 
% subplot(2,1,1); hold on;
% plot(t,x(:,1)*1e3,'linewidth',1.5); %ylabel('$l - l_0$')
% subplot(2,1,2); hold on;
% plot(x(:,2),x(:,3),'linewidth',1.5); %ylabel('$l - l_0$')

% figure(103); hold on; plot(t,x(:,9:10))
% plot(t,0.42292 + t*0,'--k','linewidth',1);
% plot(t,0.39552 + t*0,'--k','linewidth',1);

%xlabel('time (s)')

figure;plot(t,kb./kbt); hold on;
% figure(1); hold on;
% plot(t,Beta*(180/pi) );


if ShowSimulation
disp('* setting up real-time deformation render...')
mesh = Gmodel('Pneulink.stl');
mesh = Blender(mesh,'Rotate',{'z',30});
mesh = Blender(mesh,'Scale',{'uniform',0.064});

mesh.bake().render();
mesh.ground([-5e-2,5e-2,-5e-2,5e-2,0,1e-1]);
mesh.update();

FPS = 30;
Steps = round((1/FPS)/(dt));
axis on;
box on; axis on;
background('w');

for ii = 1:Steps:length(t)
    mesh.reset();
    mesh = Blender(mesh,'Rotate',{'z',30});
    mesh = Blender(mesh,'Scale',{'uniform',0.064});
    mesh = Blender(mesh,'Curve', {'PCC+',KxDeg(ii),KyDeg(ii),1+Elong(ii)});
    axis([-5e-2,5e-2,-5e-2,5e-2,0,1e-1]);
    mesh.update();
    disp(t(ii));
    pause(); 
    t(ii);
    set(gca,'linewidth',1)
    if ii == 1, disp('* press any key to play'); pause;
    end
end
end

function [t,y] = odesolve(f,t,x0)

h = 0.5*mean(diff(t));
y = zeros(length(t),length(x0));
x = x0(:);

for ii = 1:length(t)
    K1 = f(t(ii),    x);
    K2 = f(t(ii)+h,  x + h*K1);
    K3 = f(t(ii)+h,  x + h*K2);
    K4 = f(t(ii)+2*h,x + 2*h*K3);
    
    x = x + (2*h/6)*(K1 + 2*K2 + 2*K3 + K4);
    
    y(ii,:) = x(:).';
end

end

function dx = DifferentialEq(t,x)
beta0 = .05;
g   = 9.81;
r   = 0.02;
l0  = 0.064;
de  = 0.12;
db  = 5.93e-7;
kb1 = 0.42292;
kb2 = 0.39552;
kb3 = -0.21293;
ke1 = 2234;
ke2 = 1740;
ke3 = -45.55;
m   = 0.0494;


Ixx = (m*r^2)/(2*l0); 
Iyy = (m*r^2)/(4*l0); 
Izz = (m*r^2)/(4*l0); 

Ae = -28.356;
Be = -0.0018193;
Ce = 0.015346;

dx  = zeros(10,1);

l = x(1);
kx = x(2);
ky = x(3);
dl = x(4);
dkx = x(5);
dky = x(6);

q  = [l; kx; ky];
dq = [dl; dkx; dky];

M = zeros(3,3);
C = zeros(3,3);
N = zeros(3,1);
Np = zeros(3,1);
F = zeros(3,1);

[qd,dqd,ddqd] = reference(t);
% J = Jacobian(x);
% p = Position(x);
% R = rotmatPCC(x);
% 
% Rd = eye(3);
% xd = [0.05*cos(2*t);0.05*sin(t);0.03]*smoothstep(t);
% dxd = [-0.1*sin(2*t);0.05*cos(t);0]*smoothstep(t);

% G = zeros(4,4);
% Gi = G;
% G(1:3,1:3) = R;
% G(1:3,4) = p;
% Gi(1:3,1:3) = Rd.';
% Gi(1:3,4) = -Rd.'*xd;
% 
% Phi = logmapSE3(Gi*G);
% T = tmapSE3(Phi);
% 
% K = diag([1e-6,1e-6,1e-6,1,1,1]);
% 
% R = K*K*T*Phi;

L = diag([1e-4,1e-4,1e-4]);
%e = J.'*inv(J*J.')*(R);
e = (q - qd);
%J = J(4:6,1:3);
% e = J.'*inv(J*J.' + eye(3)*3e-5)*(p - xd);
% der = J.'*inv(J*J.' + eye(3)*3e-5)*(J*dq-dxd);
der = (dq-dqd);
dqr = dqd - L*e;
ddqr = ddqd - L*der;
er = der + L*e;

M(1,1) = (2*m)/(kx^2 + ky^2) - (2*m*sin(l*(kx^2 + ky^2)^(1/2)))/(l*(kx^2 + ky^2)^(3/2));
M(1,2) = (3*kx*m*sin(l*(kx^2 + ky^2)^(1/2)))/(l*(kx^2 + ky^2)^(5/2)) - (kx*m*cos(l*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^2 - (2*kx*m)/(kx^2 + ky^2)^2;
M(1,3) = (3*ky*m*sin(l*(kx^2 + ky^2)^(1/2)))/(l*(kx^2 + ky^2)^(5/2)) - (ky*m*cos(l*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^2 - (2*ky*m)/(kx^2 + ky^2)^2;
M(2,1) = M(1,2);
M(2,2) = (3*m*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 3*kx^2*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*kx^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 18*l*m*(kx^2 + ky^2)^(3/2) + 6*kx^2*l*m*(kx^2 + ky^2)^(1/2) - 3*Iyy*ky^4*l*sin(2*l*(kx^2 + ky^2)^(1/2)) + 4*Iyy*kx^4*l^4*(kx^2 + ky^2)^(3/2) + 6*Iyy*ky^4*l^2*(kx^2 + ky^2)^(1/2) + 18*Izz*ky^2*l^2*(kx^2 + ky^2)^(3/2) + 4*kx^2*l^3*m*(kx^2 + ky^2)^(3/2) - 24*Ixx*kx^2*ky^2*l*sin(l*(kx^2 + ky^2)^(1/2)) - 3*Ixx*kx^2*ky^2*l*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*kx^2*ky^2*l*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*kx^2*ky^2*l^2*(kx^2 + ky^2)^(1/2) + 4*Ixx*kx^2*ky^2*l^4*(kx^2 + ky^2)^(3/2) - 24*Izz*ky^2*l*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 3*Izz*ky^2*l*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*kx^2*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*Ixx*kx^2*ky^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*Iyy*kx^2*ky^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2))/(12*l*(kx^2 + ky^2)^(7/2));
M(2,3) = (6*kx*ky*l*m*(kx^2 + ky^2)^(1/2) - 3*kx*ky*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*kx*ky*m*sin(l*(kx^2 + ky^2)^(1/2)) - 18*Izz*kx*ky*l^2*(kx^2 + ky^2)^(3/2) + 4*kx*ky*l^3*m*(kx^2 + ky^2)^(3/2) - 12*Ixx*kx*ky^3*l*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Ixx*kx^3*ky*l*sin(l*(kx^2 + ky^2)^(1/2)) + 3*Ixx*kx^3*ky*l*sin(2*l*(kx^2 + ky^2)^(1/2)) + 12*Iyy*kx*ky^3*l*sin(l*(kx^2 + ky^2)^(1/2)) - 12*Iyy*kx^3*ky*l*sin(l*(kx^2 + ky^2)^(1/2)) + 3*Iyy*kx*ky^3*l*sin(2*l*(kx^2 + ky^2)^(1/2)) - 6*Ixx*kx^3*ky*l^2*(kx^2 + ky^2)^(1/2) + 4*Ixx*kx*ky^3*l^4*(kx^2 + ky^2)^(3/2) - 6*Iyy*kx*ky^3*l^2*(kx^2 + ky^2)^(1/2) + 4*Iyy*kx^3*ky*l^4*(kx^2 + ky^2)^(3/2) + 24*Izz*kx*ky*l*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 3*Izz*kx*ky*l*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 12*Ixx*kx*ky^3*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 12*Ixx*kx^3*ky*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 12*Iyy*kx*ky^3*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 12*Iyy*kx^3*ky*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*kx*ky*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2))/(12*l*(kx^2 + ky^2)^(7/2));
M(3,1) = M(1,3);
M(3,2) = M(2,3);
M(3,3) = (3*m*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 3*ky^2*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*ky^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 18*l*m*(kx^2 + ky^2)^(3/2) + 6*ky^2*l*m*(kx^2 + ky^2)^(1/2) - 3*Ixx*kx^4*l*sin(2*l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*kx^4*l^2*(kx^2 + ky^2)^(1/2) + 4*Ixx*ky^4*l^4*(kx^2 + ky^2)^(3/2) + 18*Izz*kx^2*l^2*(kx^2 + ky^2)^(3/2) + 4*ky^2*l^3*m*(kx^2 + ky^2)^(3/2) + 24*Ixx*kx^2*ky^2*l*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*kx^2*ky^2*l*sin(l*(kx^2 + ky^2)^(1/2)) - 3*Iyy*kx^2*ky^2*l*sin(2*l*(kx^2 + ky^2)^(1/2)) + 6*Iyy*kx^2*ky^2*l^2*(kx^2 + ky^2)^(1/2) + 4*Iyy*kx^2*ky^2*l^4*(kx^2 + ky^2)^(3/2) - 24*Izz*kx^2*l*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 3*Izz*kx^2*l*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*ky^2*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*Ixx*kx^2*ky^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*Iyy*kx^2*ky^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2))/(12*l*(kx^2 + ky^2)^(7/2));
C(1,1) = ((l*(3*dkx*kx*m*sin(l*(kx^2 + ky^2)^(1/2)) - dl*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 3*dky*ky*m*sin(l*(kx^2 + ky^2)^(1/2))))/(kx^2 + ky^2)^(5/2) + (dl*m*sin(l*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2))/l^2 - (2*dkx*kx*m*(kx^2 + ky^2)^(1/2) + 2*dky*ky*m*(kx^2 + ky^2)^(1/2) + dkx*kx*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + dky*ky*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2))/(kx^2 + ky^2)^(5/2);
C(1,2) = -(72*dkx*l^2*m*(kx^2 + ky^2)^(3/2) + 48*dkx*kx^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dkx*ky^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 3*dkx*ky^2*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 18*dkx*l*m*(kx^2 + ky^2)^(3/2) - 360*dkx*kx^2*l^2*m*(kx^2 + ky^2)^(1/2) + 8*dkx*kx^2*l^3*m*(kx^2 + ky^2)^(3/2) - 48*dkx*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) - 24*dkx*kx^4*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*kx^4*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dky*kx*ky*m*sin(l*(kx^2 + ky^2)^(1/2)) + 3*dky*kx*ky*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 30*dkx*kx^2*l*m*(kx^2 + ky^2)^(1/2) + 72*dl*kx*l^2*m*(kx^2 + ky^2)^(3/2) - 12*dkx*l*m*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(3/2) + 48*dkx*l*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) + 12*Iyy*dkx*kx^4*l^4*(kx^2 + ky^2)^(3/2) + 288*dkx*kx^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 72*dl*kx^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 72*dkx*ky^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dkx*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dkx*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 360*dky*kx*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 336*dkx*kx^2*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) - 24*dkx*kx^2*ky^2*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*kx^2*ky^2*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 360*dky*kx*ky*l^2*m*(kx^2 + ky^2)^(1/2) + 8*dky*kx*ky*l^3*m*(kx^2 + ky^2)^(3/2) + 12*Ixx*dky*kx*ky^3*l^4*(kx^2 + ky^2)^(3/2) + 12*Iyy*dky*kx^3*ky*l^4*(kx^2 + ky^2)^(3/2) - 72*dl*kx*ky^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 12*Ixx*dky*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Ixx*dky*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Iyy*dky*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 12*Iyy*dky*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*dkx*kx^2*l*m*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 48*dkx*kx^2*l*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) - 48*dl*kx*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) + 12*Ixx*dkx*kx^2*ky^2*l^4*(kx^2 + ky^2)^(3/2) - 24*dky*kx*ky^3*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dky*kx^3*ky*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dky*kx*ky^3*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dky*kx^3*ky*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 30*dky*kx*ky*l*m*(kx^2 + ky^2)^(1/2) + 12*Iyy*dkx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) - 12*Izz*dkx*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(3/2) + 48*Izz*dkx*ky^2*l^2*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) + 336*dky*kx*ky*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) - 12*Ixx*dky*kx^3*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) - 12*Iyy*dky*kx*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 12*Ixx*dkx*kx^2*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 12*dky*kx*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 48*dky*kx*ky*l*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) + 12*Izz*dky*kx*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(3/2) - 48*Izz*dky*kx*ky*l^2*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2))/(24*l^2*(kx^2 + ky^2)^(7/2));
C(1,3) = -(72*dky*l^2*m*(kx^2 + ky^2)^(3/2) + 24*dky*kx^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 3*dky*kx^2*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*dky*ky^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 18*dky*l*m*(kx^2 + ky^2)^(3/2) - 360*dky*ky^2*l^2*m*(kx^2 + ky^2)^(1/2) + 8*dky*ky^2*l^3*m*(kx^2 + ky^2)^(3/2) - 48*dky*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) - 24*dky*ky^4*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dky*ky^4*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dkx*kx*ky*m*sin(l*(kx^2 + ky^2)^(1/2)) + 3*dkx*kx*ky*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 30*dky*ky^2*l*m*(kx^2 + ky^2)^(1/2) + 72*dl*ky*l^2*m*(kx^2 + ky^2)^(3/2) - 12*dky*l*m*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(3/2) + 48*dky*l*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) + 12*Ixx*dky*ky^4*l^4*(kx^2 + ky^2)^(3/2) - 72*dky*kx^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 288*dky*ky^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 72*dl*ky^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 360*dkx*kx*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 336*dky*ky^2*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) - 24*dky*kx^2*ky^2*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dky*kx^2*ky^2*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 360*dkx*kx*ky*l^2*m*(kx^2 + ky^2)^(1/2) + 8*dkx*kx*ky*l^3*m*(kx^2 + ky^2)^(3/2) + 12*Ixx*dkx*kx*ky^3*l^4*(kx^2 + ky^2)^(3/2) + 12*Iyy*dkx*kx^3*ky*l^4*(kx^2 + ky^2)^(3/2) - 72*dl*kx^2*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 12*Ixx*dkx*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Ixx*dkx*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Iyy*dkx*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 12*Iyy*dkx*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*dky*ky^2*l*m*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 48*dky*ky^2*l*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) - 48*dl*ky*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) + 12*Iyy*dky*kx^2*ky^2*l^4*(kx^2 + ky^2)^(3/2) - 24*dkx*kx*ky^3*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*kx^3*ky*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*kx*ky^3*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*kx^3*ky*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 30*dkx*kx*ky*l*m*(kx^2 + ky^2)^(1/2) + 12*Ixx*dky*kx^4*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) - 12*Izz*dky*kx^2*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(3/2) + 48*Izz*dky*kx^2*l^2*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2) + 336*dkx*kx*ky*l^2*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) - 12*Ixx*dkx*kx^3*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) - 12*Iyy*dkx*kx*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 12*Iyy*dky*kx^2*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 12*dkx*kx*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(1/2) + 48*dkx*kx*ky*l*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(1/2) + 12*Izz*dkx*kx*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2))^2*(kx^2 + ky^2)^(3/2) - 48*Izz*dkx*kx*ky*l^2*sin((l*(kx^2 + ky^2)^(1/2))/2)^2*(kx^2 + ky^2)^(3/2))/(24*l^2*(kx^2 + ky^2)^(7/2));
C(2,1) = (48*dkx*kx^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dkx*ky^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 3*dkx*ky^2*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 8*dkx*kx^2*l^3*m*(kx^2 + ky^2)^(3/2) - 72*dl*kx*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*dkx*kx^4*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dkx*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dky*kx*ky*m*sin(l*(kx^2 + ky^2)^(1/2)) + 3*dky*kx*ky*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*dl*kx*l^2*m*(kx^2 + ky^2)^(3/2) + 12*Iyy*dkx*kx^4*l^4*(kx^2 + ky^2)^(3/2) + 6*Iyy*dkx*ky^4*l^2*(kx^2 + ky^2)^(1/2) + 18*Izz*dkx*ky^2*l^2*(kx^2 + ky^2)^(3/2) - 24*Ixx*dkx*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dkx*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 18*Izz*dky*kx*ky*l^2*(kx^2 + ky^2)^(3/2) - 24*dkx*kx^2*ky^2*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 8*dky*kx*ky*l^3*m*(kx^2 + ky^2)^(3/2) - 24*dkx*kx^2*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*dkx*kx^2*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*dl*kx*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dl*kx*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^2 - 6*Ixx*dky*kx^3*ky*l^2*(kx^2 + ky^2)^(1/2) + 12*Ixx*dky*kx*ky^3*l^4*(kx^2 + ky^2)^(3/2) - 6*Iyy*dky*kx*ky^3*l^2*(kx^2 + ky^2)^(1/2) + 12*Iyy*dky*kx^3*ky*l^4*(kx^2 + ky^2)^(3/2) - 6*Iyy*dkx*ky^4*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*Izz*dkx*ky^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Izz*dkx*ky^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 72*dl*kx*l*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 12*Ixx*dky*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Ixx*dky*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Iyy*dky*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 12*Iyy*dky*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*dkx*kx^2*ky^2*l^2*(kx^2 + ky^2)^(1/2) + 12*Ixx*dkx*kx^2*ky^2*l^4*(kx^2 + ky^2)^(3/2) - 24*dky*kx*ky^3*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dky*kx^3*ky*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 72*dl*kx*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Ixx*dkx*kx^2*ky^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*dky*kx*ky*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*dky*kx*ky*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*Izz*dky*kx*ky*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dky*kx*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Ixx*dky*kx^3*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Iyy*dky*kx*ky^3*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2))/(24*l^2*(kx^2 + ky^2)^(7/2));
C(2,2) = -(36*dkx*kx^3*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dkx*kx^3*l^4*m*(kx^2 + ky^2)^(3/2) - 8*dl*kx^2*l^3*m*(kx^2 + ky^2)^(5/2) - 9*Iyy*dky*ky^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 72*Izz*dky*ky^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 9*Izz*dky*ky^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*dkx*kx^5*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dl*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) - 6*dl*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) - 48*dl*kx^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*dl*ky^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 3*dl*ky^2*m*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 36*Izz*dky*ky*l^3*(kx^2 + ky^2)^(5/2) + 60*dkx*kx*l^2*m*(kx^2 + ky^2)^(3/2) - 8*dkx*kx*l^4*m*(kx^2 + ky^2)^(5/2) + 72*dky*ky*l^2*m*(kx^2 + ky^2)^(3/2) - 16*Iyy*dkx*kx^3*l^5*(kx^2 + ky^2)^(5/2) + 16*Iyy*dkx*kx^5*l^5*(kx^2 + ky^2)^(3/2) - 12*Iyy*dl*kx^4*l^4*(kx^2 + ky^2)^(5/2) - 24*Iyy*dky*ky^3*l^3*(kx^2 + ky^2)^(3/2) + 36*Iyy*dky*ky^5*l^3*(kx^2 + ky^2)^(1/2) - 6*Iyy*dl*ky^4*l^2*(kx^2 + ky^2)^(3/2) + 72*Izz*dky*ky^3*l^3*(kx^2 + ky^2)^(3/2) - 18*Izz*dl*ky^2*l^2*(kx^2 + ky^2)^(5/2) - 240*dkx*kx^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 120*dky*ky^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 15*dky*ky^3*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 120*Ixx*dkx*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 15*Ixx*dkx*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dkx*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dkx*kx^5*ky^2*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 120*Ixx*dky*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 15*Ixx*dky*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^2*ky^5*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 120*Iyy*dkx*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dkx*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dkx*kx^5*ky^2*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 120*Iyy*dky*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Iyy*dky*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx^2*ky^5*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 120*Izz*dkx*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 15*Izz*dkx*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Izz*dky*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 3*Izz*dky*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*dkx*kx^3*ky^2*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dky*kx^2*ky^3*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*kx*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*dkx*kx*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dl*kx^2*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dl*kx^2*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dky*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*dky*ky*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dl*kx^4*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 12*Ixx*dkx*kx*ky^2*l^3*(kx^2 + ky^2)^(3/2) - 8*Ixx*dkx*kx*ky^2*l^5*(kx^2 + ky^2)^(5/2) - 12*Ixx*dky*kx^2*ky*l^3*(kx^2 + ky^2)^(3/2) - 8*Ixx*dky*kx^2*ky*l^5*(kx^2 + ky^2)^(5/2) + 36*Iyy*dkx*kx*ky^4*l^3*(kx^2 + ky^2)^(1/2) + 16*Iyy*dky*kx^4*ky*l^5*(kx^2 + ky^2)^(3/2) + 72*Izz*dkx*kx*ky^2*l^3*(kx^2 + ky^2)^(3/2) - 72*dkx*kx*ky^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 21*dkx*kx*ky^2*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 288*dky*kx^2*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 6*dky*kx^2*ky*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 6*Iyy*dky*ky^5*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Iyy*dl*ky^4*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Izz*dky*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dky*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Izz*dl*ky^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) - 6*Izz*dl*ky^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) + 36*dky*kx^2*ky*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dky*kx^2*ky*l^4*m*(kx^2 + ky^2)^(3/2) + 168*dkx*kx^3*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dkx*kx^3*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 48*Ixx*dkx*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*dkx*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*Ixx*dky*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*dky*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 48*Iyy*dkx*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 21*Iyy*dkx*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 48*Iyy*dky*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 120*Izz*dkx*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 15*Izz*dkx*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*Izz*dky*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 6*Izz*dky*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 36*Ixx*dkx*kx^3*ky^2*l^3*(kx^2 + ky^2)^(1/2) + 16*Ixx*dkx*kx^3*ky^2*l^5*(kx^2 + ky^2)^(3/2) + 36*Ixx*dky*kx^2*ky^3*l^3*(kx^2 + ky^2)^(1/2) + 16*Ixx*dky*kx^2*ky^3*l^5*(kx^2 + ky^2)^(3/2) - 6*Ixx*dl*kx^2*ky^2*l^2*(kx^2 + ky^2)^(3/2) - 12*Ixx*dl*kx^2*ky^2*l^4*(kx^2 + ky^2)^(5/2) + 24*dky*kx^4*ky*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 168*dky*kx^2*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dky*kx^2*ky*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*dl*kx^2*ky^2*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 168*Ixx*dkx*kx^3*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Ixx*dkx*kx^3*ky^2*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 168*Ixx*dky*kx^2*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Ixx*dky*kx^2*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Ixx*dl*kx^2*ky^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 168*Iyy*dkx*kx^3*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 168*Iyy*dky*kx^2*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 48*Ixx*dkx*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 48*Ixx*dky*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 48*Iyy*dkx*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Iyy*dkx*kx*ky^4*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 48*Iyy*dky*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Izz*dkx*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dkx*kx*ky^2*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Ixx*dl*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*Ixx*dl*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*Iyy*dl*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*Iyy*dl*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2))/(24*l^2*(kx^2 + ky^2)^(9/2));
C(2,3) = -(3*Ixx*dky*kx^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dky*kx^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dky*kx^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 9*Iyy*dkx*ky^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Izz*dky*kx^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 3*Izz*dky*kx^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 72*Izz*dkx*ky^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 9*Izz*dkx*ky^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 72*Izz*dky*kx*l^3*(kx^2 + ky^2)^(5/2) - 36*Izz*dkx*ky*l^3*(kx^2 + ky^2)^(5/2) - 84*dky*kx*l^2*m*(kx^2 + ky^2)^(3/2) - 8*dky*kx*l^4*m*(kx^2 + ky^2)^(5/2) + 72*dkx*ky*l^2*m*(kx^2 + ky^2)^(3/2) + 36*Ixx*dky*kx^3*l^3*(kx^2 + ky^2)^(3/2) - 36*Ixx*dky*kx^5*l^3*(kx^2 + ky^2)^(1/2) - 8*Iyy*dky*kx^3*l^5*(kx^2 + ky^2)^(5/2) - 24*Iyy*dkx*ky^3*l^3*(kx^2 + ky^2)^(3/2) + 36*Iyy*dkx*ky^5*l^3*(kx^2 + ky^2)^(1/2) - 72*Izz*dky*kx^3*l^3*(kx^2 + ky^2)^(3/2) + 72*Izz*dkx*ky^3*l^3*(kx^2 + ky^2)^(3/2) + 168*dky*kx^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 9*dky*kx^3*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 120*dkx*ky^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 15*dkx*ky^3*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 120*Ixx*dkx*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 15*Ixx*dkx*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dkx*kx^2*ky^5*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dkx*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 96*Ixx*dky*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 120*Iyy*dkx*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Iyy*dkx*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dkx*kx^2*ky^5*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dkx*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 96*Iyy*dky*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 3*Iyy*dky*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Izz*dkx*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 3*Izz*dkx*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 168*Izz*dky*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 21*Izz*dky*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 18*Izz*dl*kx*ky*l^2*(kx^2 + ky^2)^(5/2) + 24*dkx*kx^2*ky^3*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dky*kx^3*ky^2*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 8*dl*kx*ky*l^3*m*(kx^2 + ky^2)^(5/2) - 72*dky*kx*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dky*kx*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dkx*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*dkx*ky*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 12*Ixx*dkx*kx^2*ky*l^3*(kx^2 + ky^2)^(3/2) - 8*Ixx*dkx*kx^2*ky*l^5*(kx^2 + ky^2)^(5/2) - 24*Ixx*dky*kx*ky^2*l^5*(kx^2 + ky^2)^(5/2) + 16*Ixx*dky*kx*ky^4*l^5*(kx^2 + ky^2)^(3/2) + 6*Ixx*dl*kx^3*ky*l^2*(kx^2 + ky^2)^(3/2) - 12*Ixx*dl*kx*ky^3*l^4*(kx^2 + ky^2)^(5/2) + 16*Iyy*dkx*kx^4*ky*l^5*(kx^2 + ky^2)^(3/2) + 48*Iyy*dky*kx*ky^2*l^3*(kx^2 + ky^2)^(3/2) - 72*Iyy*dky*kx*ky^4*l^3*(kx^2 + ky^2)^(1/2) + 8*Iyy*dky*kx*ky^2*l^5*(kx^2 + ky^2)^(5/2) + 6*Iyy*dl*kx*ky^3*l^2*(kx^2 + ky^2)^(3/2) - 12*Iyy*dl*kx^3*ky*l^4*(kx^2 + ky^2)^(5/2) - 144*Izz*dky*kx*ky^2*l^3*(kx^2 + ky^2)^(3/2) - 288*dkx*kx^2*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 6*dkx*kx^2*ky*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 30*dky*kx*ky^2*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Ixx*dky*kx^5*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*Iyy*dky*kx^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Iyy*dkx*ky^5*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*Izz*dky*kx^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Izz*dky*kx^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Izz*dkx*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dkx*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 36*dkx*kx^2*ky*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dkx*kx^2*ky*l^4*m*(kx^2 + ky^2)^(3/2) + 36*dky*kx*ky^2*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dky*kx*ky^2*l^4*m*(kx^2 + ky^2)^(3/2) - 24*dl*kx*ky*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 3*dl*kx*ky*m*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 48*Ixx*dkx*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*dkx*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 48*Ixx*dky*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx*ky^6*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 48*Iyy*dkx*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 48*Iyy*dky*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 18*Iyy*dky*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx*ky^6*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 48*Izz*dkx*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 6*Izz*dkx*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 144*Izz*dky*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 18*Izz*dky*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 36*Ixx*dkx*kx^2*ky^3*l^3*(kx^2 + ky^2)^(1/2) + 16*Ixx*dkx*kx^2*ky^3*l^5*(kx^2 + ky^2)^(3/2) - 72*Ixx*dky*kx^3*ky^2*l^3*(kx^2 + ky^2)^(1/2) - 36*Iyy*dky*kx^3*ky^2*l^3*(kx^2 + ky^2)^(1/2) + 16*Iyy*dky*kx^3*ky^2*l^5*(kx^2 + ky^2)^(3/2) + 24*dkx*kx^4*ky*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dky*kx*ky^4*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 168*dkx*kx^2*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dkx*kx^2*ky*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 168*dky*kx*ky^2*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dky*kx*ky^2*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 168*Ixx*dkx*kx^2*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Ixx*dkx*kx^2*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 12*Ixx*dky*kx^3*ky^2*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 168*Iyy*dkx*kx^2*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*Iyy*dky*kx^3*ky^2*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*dl*kx*ky*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dl*kx*ky*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 24*Izz*dl*kx*ky*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) + 6*Izz*dl*kx*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) + 12*Ixx*dl*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 12*Ixx*dl*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 12*Iyy*dl*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 12*Iyy*dl*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*dl*kx*ky^3*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*dl*kx^3*ky*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 48*Ixx*dkx*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 120*Ixx*dky*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 168*Ixx*dky*kx*ky^4*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*Ixx*dl*kx^3*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 48*Iyy*dkx*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 120*Iyy*dky*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 168*Iyy*dky*kx*ky^4*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 12*Iyy*dky*kx*ky^4*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*Iyy*dl*kx*ky^3*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 48*Izz*dky*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 12*Izz*dky*kx*ky^2*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2))/(24*l^2*(kx^2 + ky^2)^(9/2));
C(3,1) = (24*dky*kx^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 3*dky*kx^2*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*dky*ky^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 8*dky*ky^2*l^3*m*(kx^2 + ky^2)^(3/2) - 72*dl*ky*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*dky*ky^4*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dky*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dky*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dkx*kx*ky*m*sin(l*(kx^2 + ky^2)^(1/2)) + 3*dkx*kx*ky*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*dl*ky*l^2*m*(kx^2 + ky^2)^(3/2) + 6*Ixx*dky*kx^4*l^2*(kx^2 + ky^2)^(1/2) + 12*Ixx*dky*ky^4*l^4*(kx^2 + ky^2)^(3/2) + 18*Izz*dky*kx^2*l^2*(kx^2 + ky^2)^(3/2) + 24*Ixx*dky*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Ixx*dky*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dky*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 18*Izz*dkx*kx*ky*l^2*(kx^2 + ky^2)^(3/2) - 24*dky*kx^2*ky^2*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 8*dkx*kx*ky*l^3*m*(kx^2 + ky^2)^(3/2) - 24*dky*ky^2*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*dky*ky^2*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*dl*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dl*ky*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^2 - 6*Ixx*dkx*kx^3*ky*l^2*(kx^2 + ky^2)^(1/2) + 12*Ixx*dkx*kx*ky^3*l^4*(kx^2 + ky^2)^(3/2) - 6*Iyy*dkx*kx*ky^3*l^2*(kx^2 + ky^2)^(1/2) + 12*Iyy*dkx*kx^3*ky*l^4*(kx^2 + ky^2)^(3/2) - 6*Ixx*dky*kx^4*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*Izz*dky*kx^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Izz*dky*kx^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 72*dl*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 12*Ixx*dkx*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Ixx*dkx*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Iyy*dkx*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2)) - 12*Iyy*dkx*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Iyy*dky*kx^2*ky^2*l^2*(kx^2 + ky^2)^(1/2) + 12*Iyy*dky*kx^2*ky^2*l^4*(kx^2 + ky^2)^(3/2) - 24*dkx*kx*ky^3*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) - 24*dkx*kx^3*ky*l^2*m*sin(l*(kx^2 + ky^2)^(1/2)) + 72*dl*ky*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Iyy*dky*kx^2*ky^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*dkx*kx*ky*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*dkx*kx*ky*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*Izz*dkx*kx*ky*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dkx*kx*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Ixx*dkx*kx^3*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Iyy*dkx*kx*ky^3*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2))/(24*l^2*(kx^2 + ky^2)^(7/2));
C(3,2) = -(24*Ixx*dkx*ky^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 9*Ixx*dky*kx^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Iyy*dkx*ky^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 3*Iyy*dkx*ky^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 72*Izz*dky*kx^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 9*Izz*dky*kx^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Izz*dkx*ky^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 3*Izz*dkx*ky^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 36*Izz*dky*kx*l^3*(kx^2 + ky^2)^(5/2) + 72*Izz*dkx*ky*l^3*(kx^2 + ky^2)^(5/2) + 72*dky*kx*l^2*m*(kx^2 + ky^2)^(3/2) - 84*dkx*ky*l^2*m*(kx^2 + ky^2)^(3/2) - 8*dkx*ky*l^4*m*(kx^2 + ky^2)^(5/2) - 24*Ixx*dky*kx^3*l^3*(kx^2 + ky^2)^(3/2) + 36*Ixx*dky*kx^5*l^3*(kx^2 + ky^2)^(1/2) - 8*Ixx*dkx*ky^3*l^5*(kx^2 + ky^2)^(5/2) + 36*Iyy*dkx*ky^3*l^3*(kx^2 + ky^2)^(3/2) - 36*Iyy*dkx*ky^5*l^3*(kx^2 + ky^2)^(1/2) + 72*Izz*dky*kx^3*l^3*(kx^2 + ky^2)^(3/2) - 72*Izz*dkx*ky^3*l^3*(kx^2 + ky^2)^(3/2) - 120*dky*kx^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 15*dky*kx^3*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 168*dkx*ky^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 9*dkx*ky^3*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 96*Ixx*dkx*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 3*Ixx*dkx*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dkx*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 120*Ixx*dky*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Ixx*dky*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dky*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dky*kx^5*ky^2*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 96*Iyy*dkx*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 120*Iyy*dky*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 15*Iyy*dky*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dky*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dky*kx^5*ky^2*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 168*Izz*dkx*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 21*Izz*dkx*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Izz*dky*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 3*Izz*dky*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 18*Izz*dl*kx*ky*l^2*(kx^2 + ky^2)^(5/2) + 24*dkx*kx^2*ky^3*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dky*kx^3*ky^2*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) - 8*dl*kx*ky*l^3*m*(kx^2 + ky^2)^(5/2) + 24*dky*kx*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*dky*kx*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 72*dkx*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dkx*ky*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 48*Ixx*dkx*kx^2*ky*l^3*(kx^2 + ky^2)^(3/2) - 72*Ixx*dkx*kx^4*ky*l^3*(kx^2 + ky^2)^(1/2) + 8*Ixx*dkx*kx^2*ky*l^5*(kx^2 + ky^2)^(5/2) + 16*Ixx*dky*kx*ky^4*l^5*(kx^2 + ky^2)^(3/2) + 6*Ixx*dl*kx^3*ky*l^2*(kx^2 + ky^2)^(3/2) - 12*Ixx*dl*kx*ky^3*l^4*(kx^2 + ky^2)^(5/2) - 24*Iyy*dkx*kx^2*ky*l^5*(kx^2 + ky^2)^(5/2) + 16*Iyy*dkx*kx^4*ky*l^5*(kx^2 + ky^2)^(3/2) - 12*Iyy*dky*kx*ky^2*l^3*(kx^2 + ky^2)^(3/2) - 8*Iyy*dky*kx*ky^2*l^5*(kx^2 + ky^2)^(5/2) + 6*Iyy*dl*kx*ky^3*l^2*(kx^2 + ky^2)^(3/2) - 12*Iyy*dl*kx^3*ky*l^4*(kx^2 + ky^2)^(5/2) - 144*Izz*dkx*kx^2*ky*l^3*(kx^2 + ky^2)^(3/2) - 30*dkx*kx^2*ky*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 288*dky*kx*ky^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 6*dky*kx*ky^2*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*dky*kx^5*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 24*Ixx*dkx*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Iyy*dkx*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Iyy*dkx*ky^5*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*Izz*dky*kx^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dky*kx^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 24*Izz*dkx*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Izz*dkx*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 36*dkx*kx^2*ky*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dkx*kx^2*ky*l^4*m*(kx^2 + ky^2)^(3/2) + 36*dky*kx*ky^2*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dky*kx*ky^2*l^4*m*(kx^2 + ky^2)^(3/2) - 24*dl*kx*ky*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 3*dl*kx*ky*m*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 48*Ixx*dkx*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 18*Ixx*dkx*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dkx*kx^6*ky*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 48*Ixx*dky*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 48*Iyy*dkx*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^6*ky*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 48*Iyy*dky*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Iyy*dky*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 144*Izz*dkx*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 18*Izz*dkx*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*Izz*dky*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 6*Izz*dky*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 36*Ixx*dkx*kx^2*ky^3*l^3*(kx^2 + ky^2)^(1/2) + 16*Ixx*dkx*kx^2*ky^3*l^5*(kx^2 + ky^2)^(3/2) - 72*Iyy*dkx*kx^2*ky^3*l^3*(kx^2 + ky^2)^(1/2) + 36*Iyy*dky*kx^3*ky^2*l^3*(kx^2 + ky^2)^(1/2) + 16*Iyy*dky*kx^3*ky^2*l^5*(kx^2 + ky^2)^(3/2) + 24*dkx*kx^4*ky*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dky*kx*ky^4*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 168*dkx*kx^2*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dkx*kx^2*ky*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 168*dky*kx*ky^2*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dky*kx*ky^2*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 6*Ixx*dkx*kx^2*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 168*Ixx*dky*kx^3*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 12*Iyy*dkx*kx^2*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 168*Iyy*dky*kx^3*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Iyy*dky*kx^3*ky^2*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*dl*kx*ky*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dl*kx*ky*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 24*Izz*dl*kx*ky*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) + 6*Izz*dl*kx*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) + 12*Ixx*dl*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 12*Ixx*dl*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 12*Iyy*dl*kx*ky^5*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 12*Iyy*dl*kx^5*ky*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*dl*kx*ky^3*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*dl*kx^3*ky*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 120*Ixx*dkx*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 168*Ixx*dkx*kx^4*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 12*Ixx*dkx*kx^4*ky*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 48*Ixx*dky*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Ixx*dl*kx^3*ky*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 120*Iyy*dkx*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 168*Iyy*dkx*kx^4*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 48*Iyy*dky*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Iyy*dl*kx*ky^3*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 48*Izz*dkx*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 12*Izz*dkx*kx^2*ky*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2))/(24*l^2*(kx^2 + ky^2)^(9/2));
C(3,3) = -(36*dky*ky^3*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dky*ky^3*l^4*m*(kx^2 + ky^2)^(3/2) - 8*dl*ky^2*l^3*m*(kx^2 + ky^2)^(5/2) - 9*Ixx*dkx*kx^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 72*Izz*dkx*kx^5*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 9*Izz*dkx*kx^5*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*dky*ky^5*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dl*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) - 6*dl*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) - 24*dl*kx^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 3*dl*kx^2*m*sin(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 48*dl*ky^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 36*Izz*dkx*kx*l^3*(kx^2 + ky^2)^(5/2) + 72*dkx*kx*l^2*m*(kx^2 + ky^2)^(3/2) + 60*dky*ky*l^2*m*(kx^2 + ky^2)^(3/2) - 8*dky*ky*l^4*m*(kx^2 + ky^2)^(5/2) - 24*Ixx*dkx*kx^3*l^3*(kx^2 + ky^2)^(3/2) + 36*Ixx*dkx*kx^5*l^3*(kx^2 + ky^2)^(1/2) - 6*Ixx*dl*kx^4*l^2*(kx^2 + ky^2)^(3/2) - 16*Ixx*dky*ky^3*l^5*(kx^2 + ky^2)^(5/2) + 16*Ixx*dky*ky^5*l^5*(kx^2 + ky^2)^(3/2) - 12*Ixx*dl*ky^4*l^4*(kx^2 + ky^2)^(5/2) + 72*Izz*dkx*kx^3*l^3*(kx^2 + ky^2)^(3/2) - 18*Izz*dl*kx^2*l^2*(kx^2 + ky^2)^(5/2) - 120*dkx*kx^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 15*dkx*kx^3*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 240*dky*ky^3*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 120*Ixx*dkx*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 12*Ixx*dkx*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dkx*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dkx*kx^5*ky^2*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 120*Ixx*dky*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dky*kx^2*ky^5*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Ixx*dky*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 120*Iyy*dkx*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 15*Iyy*dkx*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^3*ky^4*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dkx*kx^5*ky^2*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 120*Iyy*dky*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 15*Iyy*dky*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dky*kx^2*ky^5*l^4*sin(l*(kx^2 + ky^2)^(1/2)) + 24*Iyy*dky*kx^4*ky^3*l^4*sin(l*(kx^2 + ky^2)^(1/2)) - 24*Izz*dkx*kx^3*ky^2*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 3*Izz*dkx*kx^3*ky^2*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 120*Izz*dky*kx^2*ky^3*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 15*Izz*dky*kx^2*ky^3*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 24*dkx*kx^3*ky^2*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dky*kx^2*ky^3*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 24*dkx*kx*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*dkx*kx*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 24*dky*ky*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*dky*ky*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dl*ky^2*l*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*dl*ky^2*l*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*dl*ky^4*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 16*Ixx*dkx*kx*ky^4*l^5*(kx^2 + ky^2)^(3/2) + 36*Ixx*dky*kx^4*ky*l^3*(kx^2 + ky^2)^(1/2) - 12*Iyy*dkx*kx*ky^2*l^3*(kx^2 + ky^2)^(3/2) - 8*Iyy*dkx*kx*ky^2*l^5*(kx^2 + ky^2)^(5/2) - 12*Iyy*dky*kx^2*ky*l^3*(kx^2 + ky^2)^(3/2) - 8*Iyy*dky*kx^2*ky*l^5*(kx^2 + ky^2)^(5/2) + 72*Izz*dky*kx^2*ky*l^3*(kx^2 + ky^2)^(3/2) - 288*dkx*kx*ky^2*l*m*sin(l*(kx^2 + ky^2)^(1/2)) - 6*dkx*kx*ky^2*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) - 72*dky*kx^2*ky*l*m*sin(l*(kx^2 + ky^2)^(1/2)) + 21*dky*kx^2*ky*l*m*sin(2*l*(kx^2 + ky^2)^(1/2)) + 6*Ixx*dkx*kx^5*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Ixx*dl*kx^4*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Izz*dkx*kx^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dkx*kx^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Izz*dl*kx^2*l^2*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) - 6*Izz*dl*kx^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(5/2) + 36*dkx*kx*ky^2*l^2*m*(kx^2 + ky^2)^(1/2) + 16*dkx*kx*ky^2*l^4*m*(kx^2 + ky^2)^(3/2) + 168*dky*ky^3*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dky*ky^3*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 48*Ixx*dkx*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 48*Ixx*dky*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 21*Ixx*dky*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*Iyy*dkx*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Iyy*dkx*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*Iyy*dky*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 6*Iyy*dky*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 48*Izz*dkx*kx*ky^4*l^2*sin(l*(kx^2 + ky^2)^(1/2)) - 6*Izz*dkx*kx*ky^4*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) - 120*Izz*dky*kx^4*ky*l^2*sin(l*(kx^2 + ky^2)^(1/2)) + 15*Izz*dky*kx^4*ky*l^2*sin(2*l*(kx^2 + ky^2)^(1/2)) + 36*Iyy*dkx*kx^3*ky^2*l^3*(kx^2 + ky^2)^(1/2) + 16*Iyy*dkx*kx^3*ky^2*l^5*(kx^2 + ky^2)^(3/2) + 36*Iyy*dky*kx^2*ky^3*l^3*(kx^2 + ky^2)^(1/2) + 16*Iyy*dky*kx^2*ky^3*l^5*(kx^2 + ky^2)^(3/2) - 6*Iyy*dl*kx^2*ky^2*l^2*(kx^2 + ky^2)^(3/2) - 12*Iyy*dl*kx^2*ky^2*l^4*(kx^2 + ky^2)^(5/2) + 24*dkx*kx*ky^4*l^3*m*sin(l*(kx^2 + ky^2)^(1/2)) + 168*dkx*kx*ky^2*l^2*m*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*dkx*kx*ky^2*l^2*m*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 24*dl*kx^2*ky^2*l^2*m*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 168*Ixx*dkx*kx^3*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 168*Ixx*dky*kx^2*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 168*Iyy*dkx*kx^3*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Iyy*dkx*kx^3*ky^2*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 168*Iyy*dky*kx^2*ky^3*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Iyy*dky*kx^2*ky^3*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) + 6*Iyy*dl*kx^2*ky^2*l^2*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 48*Ixx*dkx*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 48*Ixx*dky*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 6*Ixx*dky*kx^4*ky*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2) - 48*Iyy*dkx*kx*ky^2*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 48*Iyy*dky*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) + 24*Izz*dky*kx^2*ky*l^3*cos(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 6*Izz*dky*kx^2*ky*l^3*cos(2*l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(3/2) - 24*Ixx*dl*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) - 24*Ixx*dl*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*Iyy*dl*kx^2*ky^4*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2) + 24*Iyy*dl*kx^4*ky^2*l^3*sin(l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2))/(24*l^2*(kx^2 + ky^2)^(9/2));
N(1) = (l - l0)*(ke1 - ke2 + ke2*tanh(ke3*(l - l0))^2) - (beta0 - l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2)*(kb1 - kb2 + kb2*tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2) - (2*g*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2)/(l^2*(kx^2 + ky^2)) + (g*m*sin(l*(kx^2 + ky^2)^(1/2)))/(l*(kx^2 + ky^2)^(1/2));
N(2) = (g*kx*m*sin(l*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2) - (kx*l*(beta0 - l*(kx^2 + ky^2)^(1/2))*(kb1 - kb2 + kb2*tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2))/(kx^2 + ky^2)^(1/2) - (4*g*kx*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2)/(l*(kx^2 + ky^2)^2);
N(3) = (g*ky*m*sin(l*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2) - (ky*l*(beta0 - l*(kx^2 + ky^2)^(1/2))*(kb1 - kb2 + kb2*tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2))/(kx^2 + ky^2)^(1/2) - (4*g*ky*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2)/(l*(kx^2 + ky^2)^2);
F(1) = de*dl;
F(2) = db*dkx;
F(3) = db*dky;

ke1 = x(9);
ke2 = x(10);
kb1 = x(11);
kb2 = x(12);

Np(1) = (l - l0)*(ke1 - ke2 + ke2*tanh(ke3*(l - l0))^2) - (beta0 - l*(kx^2 + ky^2)^(1/2))*(kx^2 + ky^2)^(1/2)*(kb1 - kb2 + kb2*tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2) - (2*g*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2)/(l^2*(kx^2 + ky^2)) + (g*m*sin(l*(kx^2 + ky^2)^(1/2)))/(l*(kx^2 + ky^2)^(1/2));
Np(2) = (g*kx*m*sin(l*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2) - (kx*l*(beta0 - l*(kx^2 + ky^2)^(1/2))*(kb1 - kb2 + kb2*tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2))/(kx^2 + ky^2)^(1/2) - (4*g*kx*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2)/(l*(kx^2 + ky^2)^2);
Np(3) = (g*ky*m*sin(l*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2) - (ky*l*(beta0 - l*(kx^2 + ky^2)^(1/2))*(kb1 - kb2 + kb2*tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2))/(kx^2 + ky^2)^(1/2) - (4*g*ky*m*sin((l*(kx^2 + ky^2)^(1/2))/2)^2)/(l*(kx^2 + ky^2)^2);

tau_ = Ce*[0;x(7);x(8)];

p = zeros(3);
Aef = 0.035;
Aref = 4e-5;
dn = 2*pi/3;
H = [Aef,Aef,Aef; 
     Aref*cos(pi/2), Aref*cos(pi/2 - dn), Aref*cos(pi/2 - 2*dn);
     Aref*sin(pi/2), Aref*sin(pi/2 - dn), Aref*sin(pi/2 - 2*dn)];

om = 1/2;
 
p(1) = 90*tria(om*(t));
p(2) = 90*tria(om*(t - 2));
p(3) = 90*tria(om*(t - 4));
 
%om = 1/2;
P = 100;
% 
p(1) = P*tria(om*t);
p(2) = P*tria(om*(t - 2));
p(3) = P*tria(om*(t - 4));

Kp = diag([2500,3e-3,3e-3]);
Kd = diag([0.1,0.01,0.01]);
G = diag([1e7,1e7,5e-2,5e-2]);

tau = Kp*e + Kd*er + (M*ddqr + C*dqr) + N;
%tau =  Kp*e - Kd*dq;

Y = Slotingregressor(q);
dp = -G*Y.'*er;

ddq = M\(tau_ + tau +  - C*dq - N - F);

dx(1) = dl;
dx(2) = dkx;
dx(3) = dky;
dx(4) = ddq(1);
dx(5) = ddq(2);
dx(6) = ddq(3);
dx(7) = Ae*x(7) + Be*x(5);
dx(8) = Ae*x(8) + Be*x(6);
dx(9) = dp(1);
dx(10) = dp(2);
dx(11) = dp(3);
dx(12) = dp(4);

end

function k = stiffB(x,p)
kb3 = -0.21293;
k = p(:,3) + p(:,4).*(tanh(kb3.*x).^2 - 1);
end

function k = stiffE(x,p)
kb3 = -45.55;
k = p(:,1) + p(:,2).*(tanh(kb3.*x).^2 - 1);
end

function p = EndEffector(l,kx,ky)
kappa = sqrt(kx.^2 + ky.^2);
theta = atan2(ky,kx); 

sigma = l;
p = (1./kappa).*[(1-cos(kappa.*sigma)).*cos(theta),...
    (1-cos(kappa.*sigma)).*sin(theta),...
    sin(kappa.*sigma)];

end

function J = Jacobian(x)

l = x(1);
kx = x(2);
ky = x(3);

sigma = l;

J = zeros(6,3);

J(1,1) = 0;
J(1,2) = (kx*ky*sin(sigma*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2) - (kx*ky*sigma)/(kx^2 + ky^2);
J(1,3) = - (kx^2*sin(sigma*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2) - (ky^2*sigma)/(kx^2 + ky^2);
J(2,1) = 0;
J(2,2) = (ky^2*sin(sigma*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2) + (kx^2*sigma)/(kx^2 + ky^2);
J(2,3) = (kx*ky*sigma)/(kx^2 + ky^2) - (kx*ky*sin(sigma*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2);
J(3,1) = 0;
J(3,2) = (2*ky*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
J(3,3) = -(2*kx*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
J(4,1) = -(2*kx*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
J(4,2) = (2*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
J(4,3) = 0;
J(5,1) = -(2*ky*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
J(5,2) = 0;
J(5,3) = (2*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
J(6,1) = sin(sigma*(kx^2 + ky^2)^(1/2))/(kx^2 + ky^2)^(1/2);
J(6,2) = (kx*sigma)/(kx^2 + ky^2) - (kx*sin(sigma*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2);
J(6,3) = (ky*sigma)/(kx^2 + ky^2) - (ky*sin(sigma*(kx^2 + ky^2)^(1/2)))/(kx^2 + ky^2)^(3/2);

end

function [qd,dqd,ddqd] = reference(t)

A = 40;

qd1 = 0.07 + 0.001*sin(t);
qd2 = A*sin(t);
qd3 =A*cos(t);
dqd1 = 0.001*cos(t);
dqd2 = A*cos(t);
dqd3 = -A*sin(t);
ddqd1 = -0.001*cos(t);
ddqd2 = -A*sin(t);
ddqd3 = -A*cos(t);

qd = [qd1,qd2,qd3].';
dqd = [dqd1,dqd2,dqd3].';
ddqd = [ddqd1,ddqd2,ddqd3].';

end

function p = Position(x)

p = zeros(3,1);

l = x(1);
kx = x(2);
ky = x(3);
sigma = l;

P(1,1) = (2*kx*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
P(2,1) = (2*ky*sin((sigma*(kx^2 + ky^2)^(1/2))/2)^2)/(kx^2 + ky^2);
P(3,1) = sin(sigma*(kx^2 + ky^2)^(1/2))/(kx^2 + ky^2)^(1/2);

end

function R = rotmatPCC(x)
s = x(1);
kappa_x = x(2);
kappa_y = x(3);
kappa = sqrt(kappa_x^2 + kappa_y^2);
alpha = kappa*s;

ca = cos(alpha);
sa = sin(alpha);
va = 1-ca;

kx = -(kappa_y/kappa);
ky = (kappa_x/kappa);

R = [((kx^2)*va) + ca, (kx*ky*va), (ky*sa);
     (kx*ky*va), ((ky^2)*va) + ca, -(kx*sa);
     -(ky*sa), (kx*sa),  ca];

end

function Ad = Admap(p,R)

tx = p(1);
ty = p(2);
tz = p(3);

S = [0  -tz  ty
     tz   0  -tx
    -ty  tx   0];

Ad = zeros(6,6);
Ad(1:3,1:3) = R;
Ad(4:6,4:6) = R;
Ad(4:6,1:3) = S*R;

end

function Ad = AdInvmap(p,R)

tx = p(1);
ty = p(2);
tz = p(3);

S = [0  -tz  ty
     tz   0  -tx
    -ty  tx   0];

Ad = zeros(6,6);
Ad(1:3,1:3) = R.';
Ad(4:6,4:6) = R.';
Ad(4:6,1:3) = R*S.';

end

function Y = logmapSO3(X)
theta = acos(0.5*(trace(X) - 1.0));
a = 0.5*theta/sin(theta);

if ~isnan(a), Y = a*(X - X.');
else, Y = zeros(3,3);
end

end

function Y = tmapSO3(x)
tx = x(1);
ty = x(2);
tz = x(3);

X = [0  -tz  ty
     tz   0  -tx
    -ty  tx   0];

t = norm(x);

if (abs(t)) > 1e-6,
    a = sin(t)/t;
    b = (1.0-cos(t))/(t^2);
end

if (abs(t)) > 1e-6,
   Y = eye(3) - b*X + (1/t^2)*(1-a)*X*X;
else
   Y = eye(3); 
end

end

function Y = tmapInvSO3(x)
tx = x(1);
ty = x(2);
tz = x(3);

X = [0  -tz  ty
     tz   0  -tx
    -ty  tx   0];

t = norm(x);

if (abs(t)) > 1e-6,
    a = sin(t)/t;
    b = 2*(1.0-cos(t))/(t^2);
end

if (abs(t)) > 1e-6,
   Y = eye(3) + 0.5*X + (1/t^2)*(1-a/b)*X*X;
else
   Y = eye(3); 
end

end

function Y = logmapSE3(G)

R = G(1:3,1:3);
p = G(1:3,4);

Y = zeros(6,1);
X = logmapSO3(R);
W = [X(3,2),X(1,3),X(2,1)].';
Y(1:3) = W;

Wh = tmapInvSO3(W);
Y(4:6) = Wh.'*p;

end

function T = touplusSE3(u,o)

tx = u(1);
ty = u(2);
tz = u(3);

A = [0  -tz  ty
     tz   0  -tx
    -ty  tx   0];

tx = o(1);
ty = o(2);
tz = o(3);

B = [0  -tz  ty
     tz   0  -tx
    -ty  tx   0];

C = A*B + B*A;

t = norm(o);
a = 1;
b = 1;

if (abs(t)) > 1e-6,
    a = sin(t)/t;
    b = 2*(1.0-cos(t))/(t^2);
end

if (abs(t)) < 1e-6,
   T = -0.5*A;
else
   T = -0.5*b*A + (1/t^2)*(1-a)*C + ...
       (1/t^2)*B.'*A*((b-a)*B + ...
       (0.5*b - (1/t^2)*(3*(1-a)))*B*B);
end


end

function T = tmapSE3(x)

To = tmapSO3(x(1:3));
Tou = touplusSE3(x(4:6),x(1:3));
T = zeros(6,6);
T(1:3,1:3) = To;
T(4:6,4:6) = To;
T(1:3,4:6) = Tou;

end

function Y = Slotingregressor(x)
Y= zeros(3,2);
l0 = 0.064;
kb3 = -0.21293;
ke3 = -45.55;

l = x(1);
kx = x(2);
ky = x(3);
Y = zeros(3,4);
% Y(1,1) = l*(kx^2 + ky^2);
% Y(1,2) = l*(tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2 - 1)*(kx^2 + ky^2);
% Y(2,1) = kx*l^2;
% Y(2,2) = kx*l^2*(tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2 - 1);
% Y(3,1) = ky*l^2;
% Y(3,2) = ky*l^2*(tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2 - 1);
Y(1,1) = l - l0;
Y(1,2) = (tanh(ke3*(l - l0))^2 - 1)*(l - l0);
Y(1,3) = l*(kx^2 + ky^2);
Y(1,4) = l*(tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2 - 1)*(kx^2 + ky^2);
Y(2,1) = 0;
Y(2,2) = 0;
Y(2,3) = kx*l^2;
Y(2,4) = kx*l^2*(tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2 - 1);
Y(3,1) = 0;
Y(3,2) = 0;
Y(3,3) = ky*l^2;
Y(3,4) = ky*l^2*(tanh(kb3*l*(kx^2 + ky^2)^(1/2))^2 - 1);
end
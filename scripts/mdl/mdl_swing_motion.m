clr;
%% 
L = 100;    % length of robot
M = 5;     % number of modes
N = 100;    % number of discrete points on curve
H = 1/200;  % timesteps
FPS = 150;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature
%%
% generate nodal space
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

%%
shp = Shapes(Y,Modes,'L0',L);

shp.E    = 2.00;     % Young's modulus in Mpa
shp.Nu   = 0.49;     % Poisson ratio
shp.Rho  = 1000e-12; % Density in kg/mm^3
shp.Zeta = 0.01;      % Damping coefficient

shp.Gvec = [-9.81e3;0;0];

shp = shp.rebuild();

%%
mdl = Model(shp,'Tstep',H,'Tsim',5);

%% controller
%mdl.tau = @(M) Controller(M);

%%
mdl.q0(1)     = 0.1;
mdl = mdl.simulate(); 
%% 
figure(100);
subplot(2,1,1);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);

subplot(2,1,2);
plot(mdl.Log.t,mdl.Log.q(:,M+1:2*M),'LineW',2);
colororder(col);

%% animation
[rig] = setupRig(M,L,Modes);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end

%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup controller
function tau = Controller(mdl)
n = numel(mdl.Log.q);
t = mdl.Log.t;

w = 7;

tau        = zeros(n,1);
tau(end)   = 5*sin(w*t);
tau(end-1) = 5*sin(w*t);
% tau(n/2+1) = 9*smoothstep(t)*cos(t);
end

%% setup rig
function [rig, gmdl] = setupRig(M,L,Modes)

gmdl = Gmodel('Arm.stl');
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
     'SSSPower',0.005,'SSSRadius',5,'SSS',true);
 
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,mateplastic);
rig.g0 = SE3(roty(-pi),zeros(3,1));

rig = rig.render();
end


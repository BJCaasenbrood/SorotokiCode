clr; 
%% 
L = 100;   % length of robot
M = 5;     % number of modes
N = 50;   % number of discrete points on curve
H = 1/60; % timesteps
FPS = 30;  % animation speed

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
shp.Zeta = 0.1;      % Damping coefficient

shp = shp.rebuild();

%%
mdl = Model(shp,'Tstep',H,'Tsim',5);

%% controller
mdl.tau = @(M) Controller(M);

%%
mdl.q0(1)   = 0;
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
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
   %Y(:,ii) = pcc(X/L,ii,M);       % piece-wise constant
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup controller
function tau = Controller(mdl)
n = numel(mdl.Log.q);
t = mdl.Log.t;

%tau        = zeros(n,1);
J = mdl.Log.EL.J;
ge = SE3(mdl.Log.Phi,mdl.Log.p);
gd = SE3(roty(pi/2*t),[50,50*cos(t),50*sin(t)]);

k1 = 0.05;
k2 = 15;
lam1 = 1;
lam2 = 1;
Kp = diag([k1,k1,k1,k2,k2,k2]);

Xi = smoothstep(t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*isomse3(Xi);

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau = tau + mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q ...
    + 0.01*mdl.Log.EL.K*mdl.Log.dq;

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


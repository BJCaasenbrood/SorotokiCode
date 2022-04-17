clr; 
%% 
L = 100;   % length of robot
M = 8;     % number of modes
N = 30;    % number of discrete points on curve
H = 1/60;  % timesteps
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

%% build model class
mdl = Model(shp,'Tstep',H,'Tsim',15);

%% controller
mdl.tau = @(M) Controller(M);

%% simulate system
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);

%% animation
[rig,sph] = setupRig(M,L,Modes);
%[sph] = setupSphere();

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    sph.reset();
    sph = Blender(sph,'SE3x',gref(mdl.Log.t(ii)));
    sph = Blender(sph,'SE3',SE3(roty(-pi),[0;0;0]));
    sph.update();
    
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

t = mdl.Log.t;

J = mdl.Log.EL.J;
ge = SE3(mdl.Log.Phi,mdl.Log.p);
gd = gref(t);

k1   = 0.01;
k2   = 12;
lam1 = 1;
lam2 = 0.1;
Kp = diag([k1,k1,k1,k2,k2,k2]);

Xi = smoothstep(t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*isomse3(Xi);

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau = tau + mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;

end

function gd = gref(t)
w = 0.75;
gd = SE3(eye(3),[60,50*cos(w*t),50*sin(w*t)]);
end

%% setup rig
function [rig, sph] = setupRig(M,L,Modes)

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

sph = Gmodel(sSphere(5));

sph.Texture = diffuse(0.925);

sph.bake.render();

end



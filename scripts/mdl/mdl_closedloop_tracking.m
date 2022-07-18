clr;
%% 
L = 120;   % length of robot
M = 8;     % number of modes
N = 100;   % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature
%% shapes
X = linspace(0,1,N)';
Y = GenerateFunctionSpace(X,N,M,L);
shp = Shapes(Y,Modes,'L0',L);

%% set material properties
shp.Material = NeoHookeanMaterial(2,0.4);
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',10);
mdl = mdl.set('ResidualNorm',1e-5);
mdl = mdl.set('MaxIteration',100);

%% controller
mdl.q0 = ones(mdl.Shapes.NDim,1)*1e-3;
mdl.tau = @(M) Controller(M);

%% simulate system
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col)

%% animation
[rig, sph] = setupRig(M,L,Modes);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    sph.reset();
    sph = Blender(sph,'SE3x',gref(mdl.Log.t(ii)));
    sph = Blender(sph,'SE3',SE3(roty(pi/2),[0;0;0]));
    sph.update();
    
    axis([-.1*L .25*L -.25*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end

function gd = gref(t)
w = 0.5*pi;
gd = SE3(eye(3),[60 + 10*cos(w*t),-30*sin(w*t),-30*cos(w*t)]);
end

%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup controller
function tau = Controller(mdl)

t = mdl.Log.t;
J = mdl.Log.EL.J;

ge = SE3(mdl.Log.Phi,mdl.Log.gam);
gd = gref(t);

lam1 = 2500;
lam2 = 1e6;
Kp = diag([1e-5,1e-5,1e-5,1,1,1]);

Xi = smoothstep(t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*wedge(Xi);

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau = tau + mdl.Log.EL.fg + mdl.Log.EL.K*mdl.Log.q;
end

%% setup rig
function [rig, sph] = setupRig(M,L,Modes)

gmdl = Gmodel('Arm.stl');
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
      'SSSPower',0.2,'SSSRadius',0.25,'SSS',true);
 
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,mateplastic);
rig.g0 = SE3(roty(pi/2),zeros(3,1));

sph = Gmodel(sSphere(5));

sph.Texture = diffuse(0.925);

sph.bake.render();

rig = rig.render();
end




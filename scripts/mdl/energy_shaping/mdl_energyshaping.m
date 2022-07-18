clr; 
%% 
L = 120;   % length of robot
M = 8;     % number of modes
N = 50;   % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature

%% building model
% generate nodal space
X   = linspace(0,1,N)';
Y   = GenerateFunctionSpace(X,N,M);
shp = Shapes(Y,Modes,'L0',L);
shp.Gvec = [0,0,-9810].';

% set material properties
shp.Material = NeoHookeanMaterial(5,0.4);
shp = shp.rebuild(); 

% build model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',15,'ShowProcess',0);

I = eye(M);
mdl.InputMap = @(x) I(:,[1,2]);

%% controller
mdl.tau = @(M) Controller(M);

%% simulate system
mdl.q0(1) = -0.002;
mdl = mdl.simulate(); 

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

function gd = gref(~)
  gd = SE3(eye(3),[0,0,-40]);
end


%% setup controller
function [tau] = Controller(mdl)
%tau = zeros(size(mdl.Log.EL.G,2),1);
t = mdl.Log.t;
q = mdl.Log.q;
G = mdl.Log.EL.G;

% % T-map
[T] = CoordinateTransform(mdl);

q = mdl.Log.q;
p = mdl.Log.p;

A = T*G;

qu = annihil(A)*q;
pu = annihil(A)*p;
qa = A.'*q;
pa = A.'*p;

%qud = annihil(A)*qd;
%qad = A.'*qd;

%qdmod = T.'*[qad; qud];
qmod   = T.'*[qa; qu];

% FK kinematics on qmod
[gg,JJ] = mdl.Shapes.string(qmod);
mdlmod  = mdl.computeEL(qmod);

Tj = [(G.'*G)\(G.'); annihil(G)];
J  = JJ(:,:,end);
ge = gg(:,:,end); 

% wrench
[Fu] = ControllerWrench(t,ge);

% controller
lam1 = 450.0;
lam2 = 1e6;

dVedq = mdlmod.Log.EL.K*q + mdlmod.Log.EL.fg;
dFedq = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau   = pinv(G)*(dFedq + dVedq);

PDEcon = (annihil(G)*(dVedq));
c = PDEcon.'*PDEcon
end

function [T] = CoordinateTransform(mdl)
G = mdl.Log.EL.G;
T = [(G.'*G)\(G.'); annihil(G)];
end

function [Fu] = ControllerWrench(t,ge)
Kp = diag([1e-5,1e-5,1e-5,1,1,1]);
gd = gref(t);
Xi = smoothstep(5*t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*wedge(Xi);
end

% setup rig
function [rig, sph] = setupRig(M,L,Modes)
gmdl = Gmodel('Arm.stl','ShowProcess',0);
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
      'SSSPower',0.2,'SSSRadius',0.25,'SSS',true);
 
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,grey);
rig.g0 = SE3(roty(pi/2),zeros(3,1));

sph = Gmodel(sSphere(5));
sph.Texture = diffuse(0.925);
sph.bake.render();

rig = rig.render();
end

% functionals
function Y = GenerateFunctionSpace(X,N,M)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end




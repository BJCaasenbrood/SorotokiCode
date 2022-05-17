clr; 
%% 
L = 100;   % length of robot
M = 5;     % number of modes
N = 100;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature

%% shapes
% generate nodal space
X = linspace(0,1,N)';
Y = GenerateFunctionSpace(X,N,M);
shp = Shapes(Y,Modes,'L0',L);

%% set material properties
shp.Material = Ecoflex0050();
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',5);

%% controller
mdl.tau  = @(M) Controller(M);
mdl.Flow = @(M) IntegrateError(M);
mdl.z0   = 0;

%% simulate system
mdl.q0 = 1e-3*ones(shp.NDim,1);
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
gd = SE3(eye(3),[40,0,0]);
end

%%
function Y = GenerateFunctionSpace(X,N,M)
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

KI = 0.01;
lam1 = 25 + KI*mdl.Log.AUX.z;
lam2 = 1e4;
Kp = diag([1e-2,1e-2,1e-2,1,1,1]);

Xi = smoothstep(2*t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*wedge(Xi);

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau = tau + mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
end

function de = IntegrateError(mdl)
t = mdl.Log.t;
ge = SE3(mdl.Log.Phi,mdl.Log.gam);
gd = gref(t);
Xi = wedge(logmapSE3(ge\gd));
de = (Xi.')*Xi;
end

%% setup rig
function [rig, sph] = setupRig(M,L,Modes)
gmdl = Gmodel('Arm.stl','ShowProcess',0);

N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,base);
rig.g0 = SE3(roty(pi/2),zeros(3,1));

rig = rig.render();

sph = Gmodel(sSphere(5));
sph.Texture = diffuse(0.925);
sph.bake.render();

end




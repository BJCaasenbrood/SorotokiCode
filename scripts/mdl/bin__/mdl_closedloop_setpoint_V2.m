clr; 
%% 
L = 120;   % length of robot
M = 3;     % number of modes
N = 30;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature

%% shapes
% generate nodal space
Y   = pccspace(N,M);
shp = Shapes(Y,Modes,'L0',L);

%% set material properties
shp.Material = NeoHookeanMaterial(2,0.4);
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',15);

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
    
    axis([-.1*L L -.25*L .5*L -0.2*L 0.2*L]);
    view(30,30);
    drawnow();
end

function gd = gref(t)
gd = SE3(eye(3),[50 + 20*sin(1.25*t),0,-10]);
end

%% setup controller
function tau = Controller(mdl)
t = mdl.Log.t;
J = mdl.Log.EL.J;

[gg,JJ] = mdl.Shapes.string(mdl.Log.q);
% 
J = JJ(:,:,end);
ge = gg(:,:,end); %SE3(mdl.Log.Phi,mdl.Log.gam);
gd = gref(t);

KI = 0.03;
lam1 = 4e4 + KI*mdl.Log.AUX.z;
lam2 = 1e5;
Kp = diag([1e-5,1e-5,1e-5,1,1,1]);

Xi = smoothstep(2*t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*wedge(Xi);

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau = tau + mdl.Log.EL.fg + mdl.Log.EL.K*mdl.Log.q;
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
g1 = Gmodel('PneuLinkOne.stl','ShowProcess',0);
g2 = Gmodel('PneuLinkOne.stl','ShowProcess',0);
g3 = Gmodel('PneuLinkOne.stl','ShowProcess',0);

N = 300;
Y = pccspace(N,M);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L,'XTangent',true);

rig = rig.add(g1,g2,g3);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,0.34);
rig = rig.parent(2,0,0.33);
rig = rig.parent(2,1,0.67);
rig = rig.parent(3,0,0.66);
rig = rig.parent(3,1,1.01);

rig    = rig.texture(1,bluebase);
rig    = rig.texture(2,bluebase);
rig    = rig.texture(3,bluebase);

rig.g0 = SE3(roty(pi/2),zeros(3,1));

sph = Gmodel(sSphere(5));
sph.Texture = diffuse(0.925);
sph.bake.render();

rig = rig.render();
end




clr;
%% 
L = 360;   % length of robot
M = 2;     % number of modes
N = 100;   % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature

%% generate shapes
Y = pccspace(N,M);

shp = Shapes(Y,Modes,'L0',L);
shp.Material = NeoHookeanMaterial(25,0.4);
shp.Material.Rho  = 500e-12;
shp.Material.Zeta = 0.15;
shp.Gravity = [9.81e3; 0; 0];

shp = shp.rebuild();

%% model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',10);
mdl.Controller = @(M) Controller(M);
mdl = mdl.simulate(); 

%% animation
[rig, sph] = setupRig(M,L,Modes);

t = mdl.Log.t;

for ii = 1:fps(t,FPS):length(t)

    rig = rig.computeFK(mdl.Log.x(ii,1:2*M));
    rig = rig.update();
    
    sph.reset();
    sph = Blender(sph,'SE3x',gref(mdl.Log.t(ii)));
    sph = Blender(sph,'SE3',SE3(roty(-pi),[0;0;0]));
    sph.update();
    
    %axis([-.5*L .5*L -.5*L .5*L -1.25*L 0.1*L]);
     Bd = [ -70.0653  240.4748 -128.6664   64.5004 -350    1.0005];
     axis(Bd*1.25);
    view(140,30);
    drawnow();

end
%% setup controller
function tau = Controller(mdl)

J   = mdl.Systems{1}.Log.FK.J;
gam = mdl.Systems{1}.Log.FK.gam;
Phi = mdl.Systems{1}.Log.FK.Phi;
q   = mdl.Systems{1}.Log.q;

ge = SE3(Phi,gam);
gd = gref(mdl.t);

lam1 = 1e5;
lam2 = 1e6;
Kp = diag([0,0,0,1,1,1]);

Xi = logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*wedge(Xi);

fg = mdl.Systems{1}.Log.EL.fg;
K = mdl.Systems{1}.Log.EL.K;

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau = tau + fg + K*q;
end

function gd = gref(~)
 gd = SE3(eye(3),[125,0,-185]);
end

%% setup rig
function [rig,sph] = setupRig(M,L,Modes)

gmdl1 = Gmodel('Pneulink.stl','ShowProcess',0);
gmdl2 = Gmodel('Pneulink.stl','ShowProcess',0);
gmdl3 = Gmodel('SoftGripperRedux.stl','ShowProcess',0);

gmdl1 = Blender(gmdl1,'Loft',[1.00,0.92]);
gmdl2 = Blender(gmdl2,'Loft',[0.92,0.82]);
gmdl3 = Blender(gmdl3,'Scale',0.82);

gmdl1 = gmdl1.bake.render(); 
gmdl2 = gmdl2.bake.render(); 
gmdl3 = gmdl3.bake.render(); 
 
N = 200;
Y = pccspace(N,M);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L,'XTangent',true);

rig = rig.add(gmdl1);
rig = rig.add(gmdl2);
rig = rig.add(gmdl3);

rig = rig.parent(1,0,0);
rig = rig.parent(1,1,0.45);
rig = rig.parent(2,0,0.45);
rig = rig.parent(2,1,.9);
rig = rig.parent(3,1,.9);

rig.g0 = SE3(roty(-pi),zeros(3,1));
rig = rig.render();

sph = Gmodel(sSphere(15));
sph.Texture = diffuse(0.925);
sph.bake.render();

end


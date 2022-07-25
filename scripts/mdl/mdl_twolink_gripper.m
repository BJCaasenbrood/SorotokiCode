clr;
%% 
L = 360;   % length of robot
M = 2;     % number of modes
N = 100;   % number of discrete points on curve
H = 1/60; % timesteps
FPS = 30;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature

%% generate shapes
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

shp = Shapes(Y,Modes,'L0',L);
shp.Material = NeoHookeanMaterial(25,0.4);
shp.Material.Rho  = 500e-12;
shp.Material.Zeta = 0.1;
shp.Gvec = [9.81e3; 0; 0];

shp = shp.rebuild();

%% model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',10);
mdl = mdl.set('ResidualNorm',1e-5);
mdl = mdl.set('MaxIteration',100);

%% control law
mdl.tau  = @(M) Controller(M);

%%
mdl.q0 = 1e-2*ones(shp.NDim,1);
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:2*M),'LineW',2);
colororder(col);

%% animation
[rig, sph] = setupRig(M,L,Modes);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    sph.reset();
    sph = Blender(sph,'SE3x',gref(mdl.Log.t(ii)));
    sph = Blender(sph,'SE3',SE3(roty(-pi),[0;0;0]));
    sph.update();
    
    axis([-.5*L .5*L -.5*L .5*L -1.25*L 0.1*L]);
    view(30,30);
    drawnow();
end

%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = pcc(X/L,ii,M); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup controller
function tau = Controller(mdl)

t = mdl.Log.t;

%[gg, JJ] = mdl.Shapes.string(mdl.Log.q);
% 
J  = mdl.Log.EL.J;
ge = SE3(mdl.Log.Phi,mdl.Log.gam);%gg(:,:,end);
gd = gref(t);

lam1 = 1e-6;
lam2 = 1e-3;
Kp = diag([0,0,0,1,1,1]);

Xi = smoothstep(t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*wedge(Xi);

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
%tau = lam1*J.'*Fu;
tau = tau + mdl.Log.EL.fg + mdl.Log.EL.K*mdl.Log.q;
tau = tau*0;
end

function gd = gref(~)
 gd = SE3(eye(3),[125,0,-175]);
end

%% setup rig
function [rig,sph] = setupRig(M,L,Modes)

gmdl1 = Gmodel('Pneulink.stl','ShowProcess',0);
gmdl2 = Gmodel('Pneulink.stl','ShowProcess',0);
gmdl3 = Gmodel('SoftGripperRedux.stl','ShowProcess',0);

gmdl1 = gmdl1.bake.render(); 
gmdl2 = gmdl2.bake.render(); 
gmdl3 = gmdl3.bake.render(); 
 
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

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

rig = rig.texture(1,egg);
rig = rig.texture(2,egg);
rig = rig.texture(3,egg);
rig.g0 = SE3(roty(-pi),zeros(3,1));
rig = rig.render();

sph = Gmodel(sSphere(15));
sph.Texture = diffuse(0.925);
sph.bake.render();

end


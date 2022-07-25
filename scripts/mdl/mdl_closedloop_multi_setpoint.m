clr; 
%% 
L = 120;   % length of robot
M = 8;     % number of modes
N = 100;   % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 12;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature

%% shapes
% generate nodal space
X = linspace(0,1,N)';
Y = GenerateFunctionSpace(X,N,M);
shp = Shapes(Y,Modes,'L0',L);

%% set material properties
shp.Material = NeoHookeanMaterial(2,0.1);
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',8);

%% controller
mdl.tau  = @(M) Controller(M);

%% simulate system
mdl.q0 = 1e-4*ones(shp.NDim,1);
mdl    = mdl.simulate(); 

%% animation
h = [];
[rig, sph1, sph2] = setupRig(M,L,Modes);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    [gd1,gd2] = gref(mdl.Log.t(ii));
    
    sph1.reset();
    sph1 = Blender(sph1,'SE3x',gd1);
    sph1 = Blender(sph1,'SE3',SE3(roty(pi/2),[0;0;0]));
    sph1.update();
    
    sph2.reset();
    sph2 = Blender(sph2,'SE3x',gd2);
    sph2 = Blender(sph2,'SE3',SE3(roty(pi/2),[0;0;0]));
    sph2.update();
    
    axis([-.1*L 1.1*L -.25*L .5*L -0.45*L 0.1*L]);
    delete(h);
    h = shadowplot(5);
    
    view(30,30);
    
    drawnow();
end

function [g1,g2] = gref(~)
    g1 = SE3(rotx(-pi/6)*rotz(pi/5),[30,20,-10]);
    g2 = SE3(eye(3),[70,-20,-30]);
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
%J = mdl.Log.EL.J;

[g, J] = mdl.Shapes.string(mdl.Log.q);

ge(:,:,1) = g(:,:,40);
ge(:,:,2) = g(:,:,end);

Jb(:,:,1) = J(:,:,40);
Jb(:,:,2) = J(:,:,end);

[gd(:,:,1), gd(:,:,2)] = gref(t);

% extract positions
r = reshape(g(1:3,4,:),3,[]).' - gd(1:3,4,1).';
e = diag(r*r.');

[~,SigmaOptimal] = min(e);

% overwrite
ge(:,:,1) = g(:,:,SigmaOptimal);
Jb(:,:,1) = J(:,:,SigmaOptimal);

lam1 = 8e2;
lam2 = 1e5;
Kp = diag([1e-1,1e-1,1e-1,1,1,1]);

tau = 0;
for ii = 1:2
    Xi = smoothstep(2*t)*logmapSE3(ge(:,:,ii)\gd(:,:,ii));
    Fu = Kp*tmapSE3(Xi)*wedge(Xi);
    
    J = Jb(:,:,ii);
    Fu = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
    tau = tau + Fu;
end

tau = tau + mdl.Log.EL.fg + mdl.Log.EL.K*mdl.Log.q;

end

%% setup rig
function [rig, sph1, sph2] = setupRig(M,L,Modes)

gmdl = Gmodel('Arm.stl','ShowProcess',0);
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
    'SSSPower',0.065,'SSSRadius',0.85,'SSS',true);

N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,rot90(mateplastic,2));
rig.g0 = SE3(roty(pi/2),zeros(3,1));

sph1 = Gmodel(sTorus(0,0,0,9,2));
sph1.Texture = diffuse(0.925);
sph1.bake.render();

sph2 = Gmodel(sSphere(3),'ShowProcess',0);
sph2.Texture = diffuse(0.925);
sph2.bake.render();

rig = rig.render();

end




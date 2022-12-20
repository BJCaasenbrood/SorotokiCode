clr; 
%% 
L = 120;   % length of robot
M = 12;     % number of modes
N = 120;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature

%% shapes
% generate nodal space
Y = chebyspace(N,M);
shp = Shapes(Y,Modes,'Length',L,...
    'xia0',[0,0,0,1,0,0],...
    'Texture',softmath);

%% set material properties
shp = shp.setRadius(10);
shp = shp.setRamp(0.75);

shp.Gravity = [0;0;-9810];
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeEnd',10);

%% controller
mdl.Controller = @(M) Controller(M);

%% simulate system
mdl = mdl.simulate(); 

%% animation
sph = Gmodel(sSphere(4));
sph.Texture = diffuse(0.925);
sph.bake.render();

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)

    shp = shp.render(mdl.Log.x(ii,1:2*M));
    
    sph.reset();
    sph = Blender(sph,'SE3',gref(mdl.Log.t(ii)*0));
    %sph = Blender(sph,'SE3',SE3(roty(pi/2),[0;0;0]));
    sph.update();
    
    axis([-.1*L L -.1*L .5*L -0.25*L 0.1*L]);
    view(30,30);
    drawnow();
end

function gd = gref(t)
gd = SE3(eye(3),[50 + 30*sin(t),10*cos(t),10]);
end

%% setup controller
function tau = Controller(mdl)
t = mdl.t;
q  = mdl.Systems{1}.Log.q;
gg = mdl.Systems{1}.Log.FK.g;
JJ = mdl.Systems{1}.Log.FK.J;

ge = gg(:,:,end);
J  = admapinv(ge)*JJ(:,:,end);
gd = gref(0);

lam1 = 5e4;
lam2 = 1e4;
Kp = diag([1e-5,1e-5,1e-5,1,1,1]);

Xi = smoothstep(4*t)*logmapSE3(ge\gd);
Fu = Kp*tmapSE3(Xi)*wedge(Xi);

tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
tau = tau + mdl.Systems{1}.Log.EL.fg + mdl.Systems{1}.Log.EL.K*q;
end



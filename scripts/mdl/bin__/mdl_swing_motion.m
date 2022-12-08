clr;
%%
Omega = 4*pi; % rad/s
Ampli = 1320;  % Nmm (N-millimeter)

%% 
L = 250;    % length of robot
M = 8;      % number of modes
N = 75;     % number of discrete points on curve
H = 1/50;   % timesteps
FPS = 50;   % animation speed

%% 
Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,0,0,0,0],'L0',L);
shp.g0 = SE3(roty(pi/2),zeros(3,1));
shp.Gravity = [0;0;-9.81e3];

%% simulate model
mdl = Model(shp,'TimeStep',H,'TimeEnd',15);

mdl.Controller = @(M) Controller(M,Omega,Ampli);
mdl = mdl.simulate(); 

%% animation
[rig] = setupRig(shp);
t = mdl.Log.t;
for ii = 1:fps(t,FPS):length(t)

    rig = rig.computeFK(mdl.Log.x(ii,1:M));
    rig = rig.update();
    
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end

%% setup controller
function tau = Controller(mdl,w,A)
    tau = zeros(mdl.NDim/2,1);
    tau(1) = A*sin(w*mdl.t);
end

%% setup rig
function [rig, gmdl] = setupRig(shp)
gmdl = Gmodel('Arm.stl');

rig = Rig(@(x) shp.string(x),'Domain',shp.L0);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,softmath);
rig.g0 = SE3(roty(-pi/2),zeros(3,1));

rig = rig.render();
end


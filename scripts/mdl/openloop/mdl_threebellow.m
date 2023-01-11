clr; 
%% 
L = 60;     % length of robot
M = 1;      % number of modes
N = 60;     % number of discrete points on curve
H = 1/60;   % timesteps
FPS = 120;  % animation speed

Modes = [0,1,1,1,0,0];  % pure-XY curvature

%% shapes
% generate nodal space
Y   = pccspace(N,1);
shp = Shapes(Y,Modes,'Length',L,'xia0',[0,0,0,1,0,0]);
shp = shp.setBase(roty(-pi/2));
shp = shp.setRadius(8);

shp = shp.rebuild();

%% muscles
Muscle = @(s,N) 1.75*[0*s; cos(2*N*pi/3 + s*0); sin(2*N*pi/3 + s*0)];

for ii = 1:3
    shp = shp.addMuscle(@(s) Muscle(s,ii));
end

shp = shp.rebuild(); 

%% build model class,
mdl = Model(shp,'TimeEnd',10);

%% controller
mdl.Controller = @(M) Controller(M);

%% simulate system
mdl = mdl.simulate(); 

%% animation
rig = Rigging(shp);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)

    rig = rig.computeFK(mdl.Log.x(ii,1:3));
    rig = rig.update();
       
    axis([-.1*L .25*L -.25*L .5*L 0 L]);
    view(30,30);
    drawnow();
end


%% setup controller
function tau = Controller(mdl)
t = mdl.t;
w = 2;
A = .65;
f = @(x,k) clamp(cos(w*t - 2*(k-1)*pi/3),0,1);

tau = A*[f(t,1); f(t,2); f(t,3)];
end

function rig = Rigging(shp)
obj = Gmodel('Pneulink.stl','ShowProcess',0);
%obj = Blender(obj,'Translate',[0,0,3]);

rig = Rig(@(x) shp.string(x),'Domain',shp.Length);
rig = rig.add(obj);
rig = rig.parent(1,0,1e-6);
rig = rig.parent(1,1,1);
rig.g0 = SE3(roty(pi/2),zeros(3,1));
rig = rig.render();
end



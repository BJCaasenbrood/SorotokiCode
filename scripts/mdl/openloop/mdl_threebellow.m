clr; 
%% 
L = 60;    % length of robot
M = 1;     % number of modes
N = 30;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 120;  % animation speed

Modes = [0,1,1,1,0,0];  % pure-XY curvature

%% shapes
% generate nodal space
Y   = pccspace(N,1);
shp = Shapes(Y,Modes,'Length',L,'xia0',[0,0,0,1,0,0]);
shp = shp.setBase(roty(-pi/2));

%% set material properties
shp = shp.setRadius(8);

%% muscles
Muscle = @(s,N) 2*[0*s; sin(2*N*pi/3 + s*0); cos(2*N*pi/3 + s*0)];

for ii = 1:3
    shp = shp.addMuscle(@(s) Muscle(s,ii + 1e-3));
end

shp = shp.rebuild(); 

%% build model class,
mdl = Model(shp,'TimeEnd',10);

%% controller
mdl.Controller = @(M) Controller(M);

%% simulate system
mdl = mdl.simulate(); 

%% animation

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)

    shp = shp.render(mdl.Log.x(ii,1:3),mdl.Log.u(ii,:));
       
    axis([-.1*L .25*L -.25*L .5*L 0 L]);
    view(30,30);
    drawnow();
end


%% setup controller
function tau = Controller(mdl)
t = mdl.t;

A = -.5;

f = @(x,k) clamp(sin(t - 2*(k-1)*pi/3),0,1);

tau = A*[f(t,1); f(t,2); f(t,3)];
end



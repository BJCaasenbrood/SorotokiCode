clr; 
%% 
L = 120;   % length of robot
M = 4;     % number of modes
N = 30;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature

%% setpoint
sdf = sSphere(30,0,40,10);
sdf.show();

%% shapes
% generate nodal space
Y   = chebyspace(N,M);
shp = Shapes(Y,Modes,'Length',L,'xia0',[0,0,0,1,0,0]);

%% set material properties
shp = shp.setRadius(7);
shp = shp.setRamp(0.6);
shp = shp.setInputMap(@(x) [1;0;0;0]);
shp = shp.addContact(sdf);

shp.Material.Zeta = 0.25;
shp.Material.Rr   = 15;
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeEnd',10,'TimeStep',1/50);

%% controller
mdl.Controller = @(M) Controller(M);

%% simulate system
mdl = mdl.simulate(); 

%% animation
figure(101);


for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)

    shp = shp.render(mdl.Log.x(ii,1:M));
      
    axis([-.1*L .25*L -.25*L .5*L -0.1*L L]);
    view(30,30);
    drawnow();
end


%% setup controller
function tau = Controller(mdl)
t = mdl.t;

tau = -1500*smoothstep(0.15*t);
end



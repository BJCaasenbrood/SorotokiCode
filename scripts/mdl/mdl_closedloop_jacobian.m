clr; 
%% 
L = 100;   % length of robot
M = 6;     % number of modes
N = 60;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [M,M,0,0,0,0];  % pure-XY curvature

%% setpoint
Xd = [10;10;30];

%% shapes
% generate nodal space
Y   = chebyspace(N,M);
shp = Shapes(Y,Modes,'Length',L,'xia0',[0,0,0,0,0,1]);

%% set material properties
shp = shp.setRadius(8);
shp = shp.setRamp(0.75);
shp = shp.addGravity();

shp.Material.Zeta = 0.05;
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeEnd',10,'TimeStep',1/50);

%% controller
mdl.Controller = @(M) Controller(M,Xd);

%% simulate system
mdl = mdl.simulate(); 

%% animation
figure(101);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)

    shp = shp.render(mdl.Log.x(ii,1:2*M));
    
    if ii == 1, hold on; plot3(Xd(1),Xd(2),Xd(3),'ko'); end
      
    axis([-.1*L .25*L -.25*L .5*L -0.1*L L]);
    view(30,30);
    drawnow();
end


%% setup controller
function tau = Controller(mdl,Xd)
t = mdl.t;

sys = mdl.Systems{1};
q  = sys.Log.q;
fe = sys.Log.EL.K*q;
fg = sys.Log.EL.fg;
J  = sys.Log.FK.J(:,:,end);
x  = sys.Log.FK.g(1:3,4,end);

kp = 1e-3;

J = J(4:6,:);
tau = fe + fg + kp*J.'*(Xd - x);

end



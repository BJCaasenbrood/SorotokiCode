clr; 
%% 
L = 120;   % length of robot
M = 6;     % number of modes
N = 30;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature

%% setpoint
qd       = zeros(1,M);
qd(1:3)  = [20,-10,-30]*kpa;

%% shapes
% generate nodal space
Y   = chebyspace(N,M);
shp = Shapes(Y,Modes,'L0',L,'xia0',[0,0,0,0,0,1]);

pd = shp.FK(qd);

%% set material properties
shp = shp.setRadius(7);
shp = shp.setRamp(0.6);

shp.Material.Zeta = 0.05;
shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeEnd',10,'TimeStep',1/50);

%% controller
mdl.Controller = @(M) Controller(M,qd);

%% simulate system
mdl = mdl.simulate(); 

%% animation
figure(101);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)

    shp = shp.render(mdl.Log.x(ii,1:M));
    
    if ii == 1, hold on; fplot(pd,'--k'); end
      
    axis([-.1*L .25*L -.25*L .5*L -0.1*L L]);
    view(30,30);
    drawnow();
end


%% setup controller
function tau = Controller(mdl,qd)
t = mdl.t;

sys = mdl.Systems{1};
q   = sys.Log.q;
dq  = sys.Log.dq;
fg  = sys.Log.EL.fg;
fe  = sys.Log.EL.K*q;

k = 1e3;

Kp = kron(k,eye(sys.NJoint));
Kd = .0*Kp;

tau = (fe + fg) + Kp*(qd(:) - q) - Kd*dq;
end



clr; 
%% 
L = 100;   % length of robot
M = 2;     % number of modes
N = 200;   % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature

%% setpoint
Xd = [10;10;80];

%% shapes
% generate nodal space
NLinks = 2;
Y = affinechebyspace(N,M,NLinks);
shp = Shapes(Y,NLinks*Modes,'L0',L,'xia0',[0,0,0,0,0,1]);

%% set material properties
shp.Material = NeoHookeanMaterial(1,0.45);
shp = shp.setRadius(5);

shp = shp.setGravity();
shp.Material.Zeta = 0.2;
shp = shp.rebuild(); 

%%
shp = shp.setInputMap(@(x) [0;0;1;0]);

%% build model class
mdl = Model(shp,'TimeEnd',10,'TimeStep',1/50,'ShowProcess',true);
mdl.X0 = zeros(shp.NDim,1);

% set steady-state
mdl.X0(1:4) = [0.03198,-0.0237,0.0025,-0.0023];

%% controller
mdl.Controller = @(M) Controller(M,Xd);

%% simulate system
mdl = mdl.simulate(); 

%% animation
figure(101);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)

    shp = shp.render(mdl.Log.x(ii,1:2*M));
    
    %if ii == 1, hold on; plot3(Xd(1),Xd(2),Xd(3),'ko'); end
      
    axis([-.1*L .25*L -.25*L .5*L -L L]);
    view(30,30);
    drawnow();
end


%% setup controller
function tau = Controller(mdl,Xd)
t = mdl.t;
tau = zeros(mdl.NIn,1);
sys = mdl.Systems{1};
q  = sys.Log.q;
dq = sys.Log.dq;
M  = sys.Log.EL.M;
C  = sys.Log.EL.C;
G  = sys.Log.EL.G;
fe = sys.Log.EL.K*q;
fg = sys.Log.EL.fg;
Mi = inv(M);
% J  = sys.Log.FK.J(:,:,end);
% R  = sys.Log.FK.g(1:3,1:3,end);
% x  = sys.Log.FK.g(1:3,4,end);

K  = sys.Log.PH.Kinetic;
Ue = sys.Log.PH.Elastic;
Ug = sys.Log.PH.Gravity;

E  = K + Ue + Ug;
Er = -3.7332; % upright gravity energy

dqa = dq(3);
qa  = q(3);

Kd = 1e1;
Kv = -0e3;
Kp = -0e3;
tau = +Kd*tanh(E - Er);
%lam = E - Er + Kd*G.'*Mi*G + 1e-6;
%tau = (Kd*G.'*Mi*(C*dq + fg + fe) - Kv*dqa - Kp*qa)/lam
%tau = 1e3*sin(t);
end



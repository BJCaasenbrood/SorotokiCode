clr;
%% 
L   = 150;       % length of robot
M   = 12;         % number of modes
N   = 50;        % number of discrete points on curve
FPS = 150;       % animation speed
Xd = [80;-25;-20];

%% 
Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,M,0,0,0],'Length',L,'xia0',[0,0,0,1,0,0]);
shp = shp.setTexture(softmath);
%shp = shp.setBase(roty(pi/2));
%shp = shp.addGravity([0;0;9810]);
shp = shp.setRadius(12);
shp = shp.setRamp(0.85);

shp.Material = NeoHookeanMaterial(0.1,0.33);
shp.Material.Zeta = .5;

%% adding muscle
Nc = 1;

Muscle = @(s,N) [0*s; cos(N*pi/2 + s*0); sin(N*pi/2 + s*0)];

for ii = 1:4
    shp = shp.addMuscle(@(s) Muscle(s,ii + 1e-3));
end

Muscle = @(s,N) [0*s;cos(N*pi/2 + 2*pi*Nc*s); sin(N*pi/2+ 2*pi*Nc*s)];

for ii = 1:4
    shp = shp.addMuscle(@(s) Muscle(s,ii-1+ 1e-3));
end

Muscle = @(s,N) [0*s; cos(N*pi/2 - 2*pi*Nc*s); sin(N*pi/2 - 2*pi*Nc*s)];

for ii = 1:4
    shp = shp.addMuscle(@(s) Muscle(s,ii-1+ 1e-3));
end

shp = shp.rebuild();
%% simulate model
mdl = Model(shp,'TimeEnd',10,...
                'TimeStep',1/120,...
                'MaxIteration',5);
% 
mdl.Controller = @(M) Controller(M,Xd);
mdl = mdl.simulate(); 

%% animation
t = mdl.Log.t;
fig(101); plotpoint(Xd); hold on;
shp = shp.set('Umin',min(mdl.Log.u(end,:)));
shp = shp.set('Umax',max(mdl.Log.u(end,:)));

for ii = 1:fps(t,FPS):length(t)

    shp = shp.render(mdl.Log.x(ii,1:shp.NJoint),mdl.Log.u(ii,:));
    
    axis([-20  20  -20   20  -20   200]);
    view(30,30);
    drawnow;
end

%% 
function tau = Controller(mdl,Xd)
t = mdl.t;

sys = mdl.Systems{1};
q  = sys.Log.q;
dq = sys.Log.dq;
fe = sys.Log.EL.K*q;
fg = sys.Log.EL.fg;
G  = sys.Log.EL.G;
J  = sys.Log.FK.J(:,:,end);
x  = sys.Log.FK.g(1:3,4,end);


kp  = 1e-4;

J = J(4:6,:);

U = fe + fg + kp*J.'*(Xd - x) - kp*J.'*J*dq;
%U  = fe + kp*J.'*(((J*J.' + lam*eye(3))\(Xd - x)));
warning off;
Gplus = (G.'*G)\(G.');

tau = Gplus*U;
end

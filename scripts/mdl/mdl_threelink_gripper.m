clr;
%% 
L = 360;   % length of robot
M = 3;     % number of modes
N = 300;   % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 20;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature

%% generate shapes
X   = linspace(0,L,N)';

Y = zeros(N,M);
for ii = 1:M
   Y(:,ii) = pcc(X/L,ii,M); % pcc
end

shp = Shapes(Y,Modes,'L0',L);

%% setting materials
shp.Material     = NeoHookeanMaterial(25,0.4);
shp.Material.Rho = 250e-12;
shp.Gvec = [9.81e3; 0; 0];

shp = shp.rebuild();

%% model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',10);

%% class build + controller
qd = zeros(6,2);
qd(1,1) = +6*mm;
qd(1,2) = -6*mm;
qd(4,:) = +3*mm;
qd(5,:) = -7*mm;
qd(6,:) = -7*mm;

mdl.tau  = @(M) Controller(M,qd);

%% simulate
mdl = mdl.simulate(); 

% plotting
figure(99); clf;
subplot(1,2,1); plot(mdl.Log.t,mdl.Log.q,'Linew',2); hold on;
plot([0 5 5 10],[qd(:,1),qd(:,1),qd(:,2),qd(:,2)],'k--','Linew',1);
grid on; title('states $q$','interpreter','latex','fontsize',21)
subplot(1,2,2); plot(mdl.Log.t,mdl.Log.p,'Linew',2);
grid on; title('momenta $p = M(q)\dot{q}$','interpreter','latex',...
    'fontsize',21);
sorocolor

%% animation
figure(101); clf;
rig = setupRig(M,L,Modes);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
  
    axis([-.1*L .1*L -.1*L 0.1*L -1.2*L 0.1*L]);
    view(0,15);
    zoom(1.2);
    drawnow();
    background();
end

%% setup controller
function tau = Controller(mdl,Qd)
t   = mdl.Log.t;
q   = mdl.Log.q;
dq  = mdl.Log.dq;

if t < 5
    qd = Qd(:,1);
else
    qd = Qd(:,2);
end

KT  = mdl.Log.EL.K;
Kd  = mdl.Log.EL.M;
tau = mdl.Log.dUdq - KT*(q - qd) - Kd*dq;
end

%% setup rig
function rig = setupRig(M,L,Modes)
gmdl1 = Gmodel('Pneulink.stl','ShowProcess',0);
gmdl2 = gmdl1.copy();
gmdl3 = gmdl1.copy();

gmdl1 = Blender(gmdl1,'Loft',[1.00,0.96]);
gmdl2 = Blender(gmdl2,'Loft',[0.96,0.92]);
gmdl3 = Blender(gmdl3,'Loft',[0.92,0.88]);

gmdl4 = Gmodel('SoftGripperRedux.stl','ShowProcess',0);
gmdl4 = Blender(gmdl4,'Scale',0.88);

N = 200;
X = linspace(0,L,N)';

Y = zeros(N,M);
for ii = 1:M
   Y(:,ii) = pcc(X/L,ii,M); % pcc
end

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L,'XTangent',true);

rig = rig.add(gmdl1,gmdl2,gmdl3,gmdl4);

rig = rig.parent(1,0,0);
rig = rig.parent(1,1,0.333);
rig = rig.parent(2,0,0.333);
rig = rig.parent(2,1,.666);
rig = rig.parent(3,0,.666);
rig = rig.parent(3,1,.999);
rig = rig.parent(4,1,.999);

rig.g0 = SE3(roty(-pi),zeros(3,1));
rig = rig.render();
end


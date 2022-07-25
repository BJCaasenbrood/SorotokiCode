clr; 
%% 
L = 120;   % length of robot
T = 15;    % simulation time
M = 2;     % number of modes
N = 40;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

%% actuated states
aa = [1];
qd = [0.0379   -0.0019].';

%% building model
shp = BuildShapeLibary(N,M,L);

% build model class
mdl = Model(shp,'TimeStep',H,'TimeEnd',T,'ShowProcess',0);

I = eye(M);
mdl.InputMap = @(x) I(:,aa);

p = shp.FK(qd(:));

%% controller
mdl.tau = @(M) Controller(M,qd);

%% simulate system
mdl = mdl.simulate(); 

%%
figure(103);
subplot(1,2,1);
plot(mdl.Log.t,mdl.Log.q,'LineW',1.5); hold on;
plot(mdl.Log.t,mdl.Log.q*0 + qd.','k:','LineW',1);
xlim([0,T]);

subplot(1,2,2);
plot(mdl.Log.t,mdl.Log.p,'LineW',1.5);
xlim([0,T]);

%% animation
[rig] = setupRig(M,L);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)
  rig = rig.computeFK(mdl.Log.q(ii,:));
  rig = rig.update();
    
  axis([-.1*L .25*L -.25*L .5*L -0.25*L 0.25*L]);
  view(30,30);
  
  if ii == 1
     hold on; plot3(p(:,1),p(:,2),-p(:,3),'k--'); 
  end
  drawnow();
end

function gd = gref(~)
  gd = SE3(eye(3),[30,0,-40]);
end

%% setup controller
function [tau] = Controller(mdl,qd)
tau = zeros(size(mdl.Log.EL.G,2),1);
G = mdl.Log.EL.G;
q = mdl.Log.q;
p = mdl.Log.p;
t = mdl.Log.t;

% % T-map
[T] = CoordinateTransform(mdl);

q = mdl.Log.q;
p = mdl.Log.p;

A = T*G;

qu = annihil(A)*q;
pu = annihil(A)*p;
qa = A.'*q;
pa = A.'*p;

qud = annihil(A)*qd;
qad = A.'*qd;

qdmod = T.'*[qad; qud];
qmod  = T.'*[qa; qud];
MD0   = mdl.computeEL(q);
MDD   = mdl.computeEL(qmod);

Kp    = 1e3;
dVddq = -MDD.Log.dUdq + Kp*(qmod - qdmod) ;
tau   = -A.'*(dVddq);

pdec = annihil(G)*(MD0.Log.dUdq - dVddq);
c = pdec.'*pdec;
end

function [T] = CoordinateTransform(mdl)
G = mdl.Log.EL.G;
T = [(G.'*G)\(G.'); annihil(G)];
end

% setup rig
function [rig] = setupRig(M,L)
gmdl = Gmodel('Arm.stl','ShowProcess',0);
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
      'SSSPower',0.2,'SSSRadius',0.25,'SSS',true);
 
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M);

shp = Shapes(Y,[0,M,0,0,0,0],'L0',L);
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,grey);
rig.g0 = SE3(roty(pi/2),zeros(3,1));

rig = rig.render();
end

% functionals
function shp = BuildShapeLibary(N,M,L)

% generate nodal space
X   = linspace(0,1,N)';
Y   = GenerateFunctionSpace(X,N,M);
shp = Shapes(Y,[0,M,0,0,0,0],'L0',L);
shp.Gvec = [0,0,-9810].';

% set material properties
shp.Material = NeoHookeanMaterial(15,0.4);
shp.Material.Zeta = 0.01;
shp = shp.rebuild(); 
end

function Y = GenerateFunctionSpace(X,N,M)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end





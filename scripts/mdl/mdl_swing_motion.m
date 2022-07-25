clr;
%%
Omega = 3*pi; % rad/s
Ampli = 320;  % Nmm (N-millimeter)

%% 
L = 100;    % length of robot
M = 4;      % number of modes
N = 65;     % number of discrete points on curve
H = 1/250;  % timesteps
FPS = 30;   % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature
%%
% generate nodal space
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

%%
shp = Shapes(Y,Modes,'L0',L);
shp.Material = NeoHookeanMaterial(1,0.33);
shp.Gvec = [9.81e3;0;0];

shp = shp.rebuild();

%% simulate model
mdl = Model(shp,'TimeStep',H,'TimeEnd',15);

mdl.tau = @(M) Controller(M,Omega,Ampli);
mdl = mdl.simulate(); 

%% 
figure(100);
for ii = 1:M
    subplot(2,2,ii);
    cplot(mdl.Log.q(:,ii),mdl.Log.dq(:,ii),mdl.Log.t.^3,...
        barney(-1),'LineW',1.5);
    box on; grid on; axis square; 
    set(gca,'LineW',1.5);
end

%% animation
[rig] = setupRig(M,L,Modes);
h = gcf;
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
    set(h,'Name',['T = ',num2str(mdl.Log.t(ii),4)]);
end

%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); 
   %Y(:,ii) = pcc(X/L,ii,M); 
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup controller
function tau = Controller(mdl,w,A)
n = numel(mdl.Log.q);
t = mdl.Log.t;
tau = zeros(n,1);
tau(1) = A*sin(w*t);
end

%% setup rig
function [rig, gmdl] = setupRig(M,L,Modes)
gmdl = Gmodel('Arm.stl');
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

shp = Shapes(Y,Modes,'L0',L);

rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,softmath);
rig.g0 = SE3(roty(-pi),zeros(3,1));

rig = rig.render();
end


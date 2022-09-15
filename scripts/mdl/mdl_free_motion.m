clr;
%% 
L = 50;    % length of robot
M = 12;    % number of modes
N = 150;    % number of discrete points on curve
H = 1/200;  % timesteps
  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature
%%
% generate nodal space
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

%%
shp = Shapes(Y,Modes,'L0',L);

shp.Material = NeoHookeanMaterial(0.03,0.3);
shp.Gvec = [9.81e3;0;0];

shp = shp.rebuild();

%%
mdl = Model(shp,'TimeStep',H,'TimeEnd',5);

%%
mdl.q0(1)    = 0.1;
mdl.dq0(1+M) = 2;

mdl = mdl.simulate(); 
%% 
figure(100);
subplot(1,2,1);
plot(mdl.Log.t,mdl.Log.q(:,1:Modes(2)),'LineW',2);

subplot(1,2,1);
plot(mdl.Log.t,mdl.Log.Psi + mdl.Log.Kin + mdl.Log.Vg,'LineW',2);

%% animation
FPS = 50;
fig(101,[9,9]); clf
[rig] = setupRig(M,L,Modes);

h = [];
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    axis([-.45*L .45*L -.45*L .45*L -1.2*L 0.01*L]);
    view(30,30);
%     delete(h);
%     h = shadowplot(5);
    background(metropolis); drawnow();
    
    if ii == 1, gif('swing.gif','frame',gcf,'nodither');
    else, gif;
    end
    
    %drawnow();
end

%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup rig
function [rig, gmdl] = setupRig(M,L,Modes)

gmdl = Gmodel('Arm.stl');
% gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
%      'SSSPower',0.005,'SSSRadius',5,'SSS',true);
 
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


clr;
%% 
L = 250;   % length of robot
M = 8;     % number of modes
N = 150;   % number of discrete points on curve
H = 1/50;  % timesteps
FPS = 150;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature
%%
% generate nodal space
X = linspace(0,L,N)';
Y = chebyspace(N,M);

%% 
shp = Shapes(Y,Modes,'L0',L);
shp.Gvec = [0; 0; 0];
shp.g0 = SE3(roty(pi/2),zeros(3,1));

shp = shp.rebuild();

%%
mdl = Model(shp,'TimeEnd',30);

%%
mdl.q0(1) = 0.05;
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)
    shp = shp.render(mdl.Log.q(ii,:));
    
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end

% 
% %%
% function Y = GenerateFunctionSpace(X,N,M,L)
% % loop over functional space
% Y = zeros(N,M);
% 
% for ii = 1:M
%    Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
% end
% 
% % ensure its orthonormal (gram–schmidt)
% Y = gsogpoly(Y,X);
% end

%% setup rig
% function [rig, gmdl] = setupRig(M,L,Modes)
% 
% T = 50;
% gmdl = Gmodel('Cylinder.stl','Texture',metal);
% gmdl = Blender(gmdl,'Scale',{'z',T});
% gmdl = Blender(gmdl,'Scale',1/T);
% %gmdl = Blender(gmdl,'Fix');
%  
% gmdl = gmdl.bake.render(); 
%  
% N = 200;
% X = linspace(0,L,N)';
% Y = GenerateFunctionSpace(X,N,M,L);
% 
% shp = Shapes(Y,Modes,'L0',L);
%  
% rig = Rig(@(x) shp.string(x),'Domain',L);
% 
% rig = rig.add(gmdl);
% rig = rig.parent(1,0,0);
% rig = rig.parent(1,1,1);
% 
% rig.g0 = SE3(roty(-pi),zeros(3,1));
% 
% rig = rig.render();
% end


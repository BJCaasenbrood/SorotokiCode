clr;
%% parameters
M   = 10;
Y   = chebyspace(60,M);

L   = 500;       % mm
Rho = 7850e-12;  % kg/mm^3
E0  = 210e3;     % MPa
Nu0 = 0.3;       % (-)
R   = 1;         % mm

shp = Shapes(Y,[0,M,M,0,0,0],'Length',L,'xia0',[0;0;0;1;0;0]);
shp.Material = NeoHookeanMaterial(E0,Nu0);

shp.Material.Zeta = 1e-6;
shp.Material.Rho  = Rho;

shp = shp.setRadius(R);
shp = shp.rebuild();

mdl = Model(shp,'TimeEnd',3,'TimeStep',1/50);
mdl.Controller = @(x) Control(x);

%%
mdl = mdl.simulate();

%%
t = mdl.Log.t;
Ndof = shp.NJoint;
P = [];
for ii = 1:numel(t)
   p = shp.FK(mdl.Log.x(ii,1:Ndof));
   P = [P;p(end,:)];   
end
    
%%
fig(103); cla;
load('matlab_grazioso.mat');

plot(X(:,1),X(:,2),'--','LineW',2,'Color',color_gray_darkest); hold on;
plot(Y(:,1),Y(:,2),'-.','LineW',2,'Color',color_gray_darkest); 
plot(Z(:,1),Z(:,2),':','LineW',2,'Color',color_gray_darkest); 
plot(t,P(:,1)/1e3,'LineW',2,'Color',col(1)); hold on;
plot(t,P(:,2)/1e3,'LineW',2,'Color',col(2)); hold on;
plot(t,P(:,3)/1e3,'LineW',2,'Color',col(3)); hold on;

%%

function tau = Control(mdl)
t = mdl.t;
J = mdl.Systems{1}.Log.FK.J(:,:,end);
g = mdl.Systems{1}.Log.FK.g(:,:,end);

% inertial frame Jacobian
Je = admapinv(g(1:3,1:3))*J;

F = 0.14*t;
M = 490*t;

W = [0;0;M;0;0;F];
tau = (Je).'*W;
end
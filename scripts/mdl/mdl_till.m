clr;
%% parameters
M   = 3;
Y   = chebyspace(100,M);

L   = 517;       % mm
Rho = 7800e-12;  % kg/mm^3
E0  = 200e3;     % MPa
Nu0 = 0.3;       % (-)
R   = 1.42;      % mm

shp = Shapes(Y,[0,M,M,0,0,0],'L0',L);
shp.Material = NeoHookeanMaterial(E0,Nu0);

shp.Material.Zeta = 1e-6;
shp.Material.Rho  = Rho;

shp = shp.setRadius(R);
shp = shp.rebuild();

mdl = Model(shp,'TimeEnd',3,'TimeStep',1/1e3);
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
fig(103);
% load('matlab_grazioso.mat');
% 
% plot(X(:,1),X(:,2),'LineW',2,'Color',color_gray_light); hold on;
% plot(Y(:,1),Y(:,2),'LineW',2,'Color',color_gray_light); 
% plot(Z(:,1),Z(:,2),'LineW',2,'Color',color_gray_light); 
% plot(t,P(:,1)/1e3,'LineW',2,'Color',col(1)); hold on;
% plot(t,P(:,2)/1e3,'LineW',2,'Color',col(2)); hold on;
plot(t,P(:,3)/1e3,'LineW',2,'Color',col(3)); hold on;

%%
% figure(101);
% shp = shp.setRadius(3*R);
% 
% for ii = 1:numel(t)
%     shp = shp.render(mdl.Log.x(ii,1:Ndof));
%     drawnow();    
%     view(30,30);
%     axis([-0,.5,-.5,.5,0,.5]*500);
% end

function tau = Control(mdl)
M = 50;
d = 0.016;
t = mdl.t;
J = mdl.Systems{1}.Log.FK.J(:,:,30);

if t < 0.5*d
   F = M*(t/(0.5*d));
elseif (t >= 0.5*d) && (t <= d)
   F = M*(2 - t/(0.5*d)); 
else
   F = 0;
end

W = [0;0;0;0;0;F];
tau = (J).'*W;
end
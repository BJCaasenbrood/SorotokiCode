clr;
%% generate environment
sdf = Environment;
sdf.show();
colormap(barney);
% obj = Gmodel(sdf,'Quality',100);
% obj.bake.render();

%% generate soft beam
shp = BuildShapeLibary(250,8,10);
nq = shp.NDim;


p = shp.FK(rand(nq,1)*1e-1);
fplot(p,'LineW',4);

%% sdf func
function sdf = Environment()
X1 = 1; X2 = 4;
Z1 = 4; Y1 = 1;
Z2 = 5; Y2 = -1;

W  = 5;
DT = 0.25;
O = sSphere(0.25);
C1 = sCube(X1,X1+DT,-W,W,0,2*W); S1 = sSphere(X1,Y1,Z1,1);
C2 = sCube(X2,X2+DT,-W,W,0,2*W); S2 = sSphere(X2,Y2,Z2,1);
sdf = (C1-S1) + (C2-S2);
end

% functionals
function shp = BuildShapeLibary(N,M,L)

% generate nodal space
X   = linspace(0,1,N)';
Y   = GenerateFunctionSpace(X,N,M);
shp = Shapes(Y,[0,M,M,0,0,0],'L0',L);
shp.Gvec = [0,0,-9810].';
shp.g0 = SE3(roty(-pi/2),[0;0;0]);

% set material properties
shp.Material = NeoHookeanMaterial(2,0.4);
shp.Material.Zeta = 0.05;
shp = shp.rebuild(); 
end

function Y = GenerateFunctionSpace(X,N,M)
Y = zeros(N,M);
for ii = 1:M
   Y(:,ii) = chebyshev(X,ii-1); % chebyshev
end
Y = gsogpoly(Y,X);
end
close all; clc; clear;
%% set signed distance function
sdf = @(x) Square(x,10,20);
%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,20],'NElem',500);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'Penal',3,...
              'FilterRadius',1e-3,...
              'VolumetricPressure',true);
              
%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,1]);
fem = fem.AddConstraint('PressureCell',fem.FindElements(...
                        'Location',[5,5],1),[4e-3,0]);

fem.Material = Ecoflex0030(1e36);

%% generate void region
f = @(x) dRectangle(x,2,8,3.75,20);
d = f(fem.Center); id = find(d(:,end) < 0);
fem.Density(id) = 1e-2;

%% solving
fem.solve();

function D = Square(P,W,H)
D = dRectangle(P,0,W,0,H);
end
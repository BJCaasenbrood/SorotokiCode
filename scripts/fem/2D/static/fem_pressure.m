close all; clc; clear;
%% set signed distance function
sdf = sRectangle(0,10,0,20);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,20],'Quads',150);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/30,'ResidualNorm',1e-3,...
              'Nonlinear',true,'Penal',1,...
              'FilterRadius',1e-3,'VolumetricPressure',true,...
              'Linestyle','-');
              
%% add constraint
fem = fem.addSupport(fem.FindNodes('Bottom'),[1,1]);
%fem = fem.addSupport(fem.FindNodes('SW'),[1,1]);
%fem = fem.addSupport(fem.FindNodes('Top'),[0,1]);
fem = fem.addMyocyte(fem.FindElements('Location',[5,5],1),[-1*kpa,0]);

fem.Material = Ecoflex0030();

%% generate void region
f = @(x) dRectangle(x,2,8,3.75,18);
d = f(fem.Center); 
id = find(d(:,end) < 0);
fem.Density(id) = 0.01;

%% solving
fem.solve();
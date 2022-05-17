close all; clc; clear;
%% set signed distance function
sdf = sRectangle(0,10,0,20);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,20],'NElem',500);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/25,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'Penal',3,...
              'FilterRadius',1e-3,...
              'VolumetricPressure',true);
              
%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,1]);
fem = fem.AddConstraint('PressureCell',fem.FindElements(...
                        'Location',[5,5],1),[2*kpa,0]);

fem.Material = Ecoflex0050(30);

%% generate void region
f = @(x) dRectangle(x,2,8,3.75,20);
d = f(fem.Center); 
id = find(d(:,end) < 0);
fem.Density(id) = 0.25;

%% solving
fem.solve();
clear; close all; clc;

%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,3);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,5,0,3],'NElem',500);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.2,...
              'Penal',4,...
              'VolumetricPressure',true,...
              'FilterRadius',1.5,...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Support',id,[0,1]);
id = fem.FindNodes('Left');
fem = fem.AddConstraint('Support',id,[1,0]);

id = fem.FindNodes('NW');
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,1]);
id = fem.FindNodes('SE');
fem = fem.AddConstraint('Output',id,[1,0]);
fem = fem.AddConstraint('Spring',id,[1,0]);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[-1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0],1);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();
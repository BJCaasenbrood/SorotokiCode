clear; close all; clc;

%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,7.5);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,7.5],...
              'NElem',400,...
              'MaxIteration',50,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.2,...
              'Penal',4,...
              'Nonlinear',false,...
              'FilterRadius',1,...
              'MaxIterationMMA',100,...
              'VolumetricPressure',true,...
              'ReflectionPlane',[1,1],...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant');

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[-1e-3,0]);

id = fem.FindNodes('NW');
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1e-3]);
% % 
% id = fem.FindNodes('SE');
% fem = fem.AddConstraint('Output',id,[-1,-1]);
% fem = fem.AddConstraint('Spring',id,[1e-3,0]);
%% material
fem.Material = Dragonskin10A;

%% set density
fem = fem.initialTopology('Hole',[0,0],3);

%% solving
fem.optimize();

%% reconstruct
fem = fem.set('ReflectionPlane',[1,1],...
        'CellRepetion',[2,6]);

fem = fem.former(10);
fem.showISO(0.15,.5);

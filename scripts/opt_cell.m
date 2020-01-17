clear; close all; clc;

%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,10],...
              'NElem',900,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.2,...
              'Penal',1,...
              'PenalMax',5,...
              'Nonlinear',false,...
              'FilterRadius',1,...
              'MaxIterationMMA',100,...
              'VolumetricPressure',true,...
              'ReflectionPlane',[1,1],...
              'Nonlinear',true,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,0]);
id = fem.FindNodes('Bottom'); 
fem = fem.AddConstraint('Support',id,[0,1]);

% id = fem.FindNodes('Location',[15,5],1); 
% fem = fem.AddConstraint('Output',id,[-1,0]);
% fem = fem.AddConstraint('Spring',id,[0.05,0]);

id = fem.FindElements('Location',[1,1],1);
fem = fem.AddConstraint('PressureCell',id,[0.01e-3,0]);

id = fem.FindNodes('Location',[1,10],1); 
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,05]);

%% material
fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',1.5,'D2',2.0,'D3',1.0);

%fem.Material = MooneyMaterial('C10',3,'K',50);

%fem.Material = LinearMaterial('E',3,'Nu',0.49);

%% set density
fem = fem.initialTopology('Hole',[0,0],5);

%% solving
fem.optimize();

%% reconstruct
fem.set('ReflectionPlane',[1,1],...
        'CellRepetion',[4,12]);

fem.former(10);
fem.showISO(0.35,.5);

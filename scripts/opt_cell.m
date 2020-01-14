clear; close all; clc;

%% set signed distance function
sdf = @(x) dRectangle(x,0,15,0,10);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,15,0,10],...
              'NElem',500,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.3,...
              'Penal',1,...
              'PenalMax',4,...
              'Nonlinear',false,...
              'FilterRadius',1,...
              'MaxIterationMMA',100,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,0]);
id = fem.FindNodes('Bottom'); 
fem = fem.AddConstraint('Support',id,[0,1]);

id = fem.FindNodes('Location',[15,10*0.2],1); 
fem = fem.AddConstraint('Load',id,[1,0]);
fem = fem.AddConstraint('Spring',id,[1,0]);

id = fem.FindNodes('Location',[15*0.33,10],1); 
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

%% material
 fem.Material = YeohMaterial('C1',1,'C2',0.1,'C3',0.1,...
     'D1',1,'D2',1,'D3',1);

%fem.Material = LinearMaterial('E',3,'Nu',0.49);

%% set density
fem = fem.initialTopology([1,1],2);

%% solving
fem.optimize();

%% reconstruct
fem = fem.former;
fem.showTopo;

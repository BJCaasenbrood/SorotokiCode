clear; close all; clc;

%% set signed distance function
sdf = @(x) PneuGrip(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[-3,13,0,5],...
              'NElem',1000,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.25,...
              'VolumetricPressure',true,...
              'FilterRadius',0.75,...
              'Nonlinear',true,...
              'MaxIterationMMA',80,...
              'ReflectionPlane',[0 -1],...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Top'); fem = fem.AddConstraint('Support',id,[0,1]);
id = fem.FindNodes('Line',[]); fem = fem.AddConstraint('Support',id,[0,1]);

id = fem.FindNodes('Location',[13,3],1);
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,10]);

id = fem.FindElements('Location',[-3,5],1);
fem = fem.AddConstraint('PressureCell',id,[0.1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[-3,5],1);

%% material
fem.Material = Ecoflex0030;

%% solving
fem.optimize();

function Dist = PneuGrip(P)
  R1 = dRectangle(P,-3,13,0,5);
  R2 = dRectangle(P,10,13,3,7);
  R3 = dRectangle(P,10,13,0,1);
  C2 = dCircle(P,10,5,2);
  C3 = dCircle(P,10,0,1);
  Dist = dDiff(dDiff(dDiff(dDiff(R1,C2),R2),R3),C3);
end

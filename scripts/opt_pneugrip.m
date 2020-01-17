clear; close all; clc;

%% set signed distance function
sdf = @(x) PneuGrip(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[-3,13,0,10],...
              'NElem',900,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.3,...
              'Penal',1,...
              'PenalMax',4,...
              'PrescribedDisplacement',false,...
              'VolumetricPressure',true,...
              'FilterRadius',0.75,...
              'Nonlinear',false,...
              'MaxIterationMMA',80,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[13,3],1);
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindNodes('Location',[13,7],1);
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindNodes('Location',[8,5],1);
fem = fem.AddConstraint('Output',id,[-1,0]);
fem = fem.AddConstraint('Spring',id,[0,0]);

id = fem.FindElements('Location',[-3,5],1);
fem = fem.AddConstraint('PressureCell',id,[-0.1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[-3,5],1);

%% material
fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
   'D1',1.5,'D2',2.0,'D3',1.0);

%fem.Material = MooneyMaterial('C10',3,'K',50);

%fem.Material = LinearMaterial('E',3,'Nu',0.49);

fem.show('E');

%% solving
fem.optimize();

%% former
fem.former(6);
fem.showSTL(0.1);


function Dist = PneuGrip(P)
  R1 = dRectangle(P,-3,13,0,10);
  R2 = dRectangle(P,10,13,3,7);
  R3 = dRectangle(P,10,13,0,1);
  R4 = dRectangle(P,10,13,9,10);
  C2 = dCircle(P,10,5,2);
  C3 = dCircle(P,10,0,1);
  C4 = dCircle(P,10,10,1);
  Dist = dDiff(dDiff(dDiff(dDiff(dDiff(dDiff(R1,C2),R2),R3),R4),C3),C4);
end


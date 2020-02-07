clear; close all; clc;

%% set signed distance function
sdf = @(x) PneuRot(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,5],...
              'NElem',700,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.2,...
              'Penal',1,...
              'PenalMax',4,...
              'PrescribedDisplacement',false,...
              'VolumetricPressure',true,...
              'FilterRadius',0.4,...
              'Nonlinear',false,...
              'ReflectionPlane',[1,1],...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Location',[0,5],1);
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Location',[5,0],1);
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[2,0],1);
fem = fem.AddConstraint('Support',id,[1,0]);
fem = fem.AddConstraint('Output',id,[-1,-1]);
fem = fem.AddConstraint('Spring',id,[1,1]);

id = fem.FindNodes('Location',[0,2],1);
fem = fem.AddConstraint('Support',id,[0,1]);
fem = fem.AddConstraint('Output',id,[1,-1]);
fem = fem.AddConstraint('Spring',id,[1,1]);
% % 
% 
% id = fem.FindElements('Location',[0,1],1);
% fem = fem.AddConstraint('Output',id,[1,0]);
% fem = fem.AddConstraint('Spring',id,[1,0]);

id = fem.FindElements('Location',[2.5,2.5],1);
fem = fem.AddConstraint('PressureCell',id,[-0.05e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[2.5,2.5],0.25);

%% material
fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',1.5,'D2',2.0,'D3',1.0);

%% solving
fem.optimize();

fem.former(3);
fem.showSTL(0.15);

function Dist = PneuRot(P)
  R1 = dRectangle(P,0,5,0,5);
  C1 = dCircle(P,0,0,5);
  C2 = dCircle(P,0,0,2);
  Dist = dIntersect(R1,dDiff(C1,C2));
end


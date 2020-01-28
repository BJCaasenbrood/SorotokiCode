clear; close all; clc;

%% set signed distance function
sdf = @(x) PneuRot(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,10],...
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
              'VolumeInfill',0.2,...
              'Penal',1,...
              'PenalMax',4,...
              'PrescribedDisplacement',false,...
              'VolumetricPressure',true,...
              'FilterRadius',0.35,...
              'Nonlinear',false,...
              'ReflectionPlane',[1,1],...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Location',[0,5],1);
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Location',[0,2],1);
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Support',id,[1,0]);
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindNodes('Left');
fem = fem.AddConstraint('Support',id,[0,1]);
fem = fem.AddConstraint('Load',id,[1,0]);
fem = fem.AddConstraint('Spring',id,[1,1]);

%% set density
fem = fem.initialTopology('Hole',[2.5,2.5],0.5);

%% material
% fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
%    'D1',1.5,'D2',2.0,'D3',1.0);
% 
fem.Material = LinearMaterial('E',3,'Nu',0.49);

fem.show('E');

%% solving
fem.optimize();

function Dist = PneuRot(P)
  R1 = dRectangle(P,0,5,0,5);
  C1 = dCircle(P,0,0,5);
  C2 = dCircle(P,0,0,2);
  Dist = dIntersect(R1,dDiff(C1,C2));
end




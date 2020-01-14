clear; close all; clc;

%% set signed distance function
sdf = @(x) PneuNet(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,20],...
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
              'VolumeInfill',0.3,...
              'Penal',1,...
              'PenalMax',4,...
              'PrescribedDisplacement',false,...
              'VolumetricPressure',true,...
              'FilterRadius',1.5,...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindNodes('Bottom'); 
fem = fem.AddConstraint('Spring',id,[1,0]);

id = fem.FindElements('Location',[10,10],1);
fem = fem.AddConstraint('PressureCell',id,[-0.01e-3,0]);

%% set density
fem = fem.initialTopology([1,1],5);

%% material
% fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
%    'D1',1.5,'D2',2.0,'D3',1.0);

fem.Material = MooneyMaterial('C10',3,'K',50);

%fem.Material = LinearMaterial('E',3,'Nu',0.49);

fem.show('E');
colormap(bluesea);
%% solving
fem.optimize();

function Dist = PneuNet(P)
  R1 = dRectangle(P,0,20,0,20);
  R2 = dRectangle(P,-4,1,4,20);
  R3 = dRectangle(P,19,24,4,20);
  C1 = dCircle(P,0,4.5,1);
  C2 = dCircle(P,20,4.5,1);
  Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end


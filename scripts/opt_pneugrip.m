clear; close all; clc;

%% set signed distance function
sdf = @(x) PneuGrip(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,13,0,10],...
              'NElem',1000,...
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
              'FilterRadius',0.75,...
              'Nonlinear',false,...
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

id = fem.FindElements('Location',[0,5],1);
fem = fem.AddConstraint('PressureCell',id,[-0.1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[0,5],1);

%% material
fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
   'D1',1.5,'D2',2.0,'D3',1.0);

%fem.Material = MooneyMaterial('C10',3,'K',50);

%fem.Material = LinearMaterial('E',3,'Nu',0.49);

fem.show('E');

%% solving
fem.optimize();

%% interpolate
verts = fem.Center; 
x = verts(:,1); 
y = verts(:,2);
BdBox = msh.BdBox;
xq = linspace(BdBox(1),BdBox(2),150); 
yq = linspace(BdBox(3),BdBox(4),150);
[xq,yq] = meshgrid(xq,yq);
P = fem.get('SpatialFilter');
vq = griddata(x,y,P*fem.Density,xq,yq);
contour(xq,yq,vq,[0.47 0.47]);

function Dist = PneuGrip(P)
  R1 = dRectangle(P,0,13,0,10);
  R2 = dRectangle(P,10,13,3,7);
  C2 = dCircle(P,10,5,2);
  Dist = dDiff(dDiff(R1,C2),R2);
end


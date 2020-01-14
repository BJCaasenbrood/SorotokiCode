clear; close all; clc;

%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,3);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,3],...
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
              'FilterRadius',0.4,...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Support',id,[0,1]);
id = fem.FindNodes('Right');
fem = fem.AddConstraint('Output',id,[-1,0]);
fem = fem.AddConstraint('Spring',id,[1,0]);

id = fem.FindNodes('Left');
fem = fem.AddConstraint('Support',id,[1,0]);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[0.05e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0],0.25);

%% material
fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',1.5,'D2',2.0,'D3',1.0);
% 
%fem.Material = LinearMaterial('E',30,'Nu',0.49);

fem.show('E');

%% solving
fem.optimize();

fem = fem.former;
fem.showTopo;

% %% reconstruct
% id = (fem.Density >= 1e-3);
% verts = fem.Center(id,:); 
% 
% Dist = @(P)dMetasphere2D(P,verts(:,1),verts(:,2),0.75,10);
% BdBox = msh.BdBox;
% x = linspace(BdBox(1),BdBox(2),100); y = linspace(BdBox(3),BdBox(4),100);
% [X,Y] = meshgrid(x,y); P = [X(:),Y(:)];
% d = Dist(P);
% 
% cla;
% D = reshape(d(:,end),[100 100]);
% hold on;
% surf(X,Y,D);
%contour(X,Y,D,[0 0])

function Dist = PneuRot(P)
  R1 = dRectangle(P,0,5,0,5);
end


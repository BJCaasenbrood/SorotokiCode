clc; clear; close all;
%% set signed distance function
L = 40;
W = 3;

sdf = @(x) dRectangle(x,0,L,0,W);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,L,0,W],...
              'NElem',150,...
              'MaxIteration',500,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% generate fem model
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/50,...
              'ResidualNorm',1e-3,...
              'Nonlinear',false,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('NW'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('NE'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[-30,0]);

fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
     'D1',15,'D2',10,'D3',10);

%% solving
fem.solve();

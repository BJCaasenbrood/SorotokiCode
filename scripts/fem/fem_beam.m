clc; clear; close all;
%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,2);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,2],...
              'NElem',500);
      
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/15,...
              'ResidualNorm',1e-3);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[1,0]);
fem.Material = Ecoflex0030;

%% solving
fem.solve();

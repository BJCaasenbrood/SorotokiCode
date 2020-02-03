clc; clear; close all;
%% set signed distance function
sdf = @(x) dRectangle(x,0,8,0,1);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,8,0,1],...
              'NElem',150);
      
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/25);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[0,-1e-3]);

%% select material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

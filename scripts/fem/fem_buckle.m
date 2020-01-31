clc; clear; close all;
%% set signed distance function
sdf = @(x) dRectangle(x,0,40,0,2);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,40,0,2],...
              'NElem',150);
      
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/50,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('NE'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[-30,0]);

fem.Material = Ecoflex0030;

%% solving
fem.solve();

clc; clear; close all;
%% set signed distance function
sdf = @(x) dRectangle(x,0,1,0,1);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,1,0,1],...
              'Center',[0.4 0.5; 0.6 0.5]);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/25,...
              'ResidualNorm',1e-9,...
              'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[7,0]);

fem.Material = Dragonskin10A;

%% solving
fem.solve();
close all; clc; clear;
%% set signed distance function
sdf = @(x) SDF(x,20,1);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,10],...
              'NElem',500);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/150,...
              'ResidualNorm',1e-3,...
              'DisplaceNorm',1e-3,...
              'Nonlinear',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[1,1]);
fem = fem.AddConstraint('Pressure',fem.FindNodes('Line',[0 20 0 0]),[1e-1,0]);
%% material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

function Dist = SDF(P,R,t)
Dist = dRectangle(P,0,R,0,t);
end

clc;  clear; close all;
%% set signed distance function
sdf = @(x) SquareHole(x,2,.5);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,5],...
              'NElem',250,...
              'MaxIteration',1e3);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/50,...
              'ResidualNorm',1e-4,...
              'DisplaceNorm',1e-4,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,1]);

%% material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

function D = SquareHole(P,H,R)
R1 = dRectangle(P,0,H/2,0,H/2);
C1 = dCircle(P,0,0,R);
D = dDiff(R1,C1);
end
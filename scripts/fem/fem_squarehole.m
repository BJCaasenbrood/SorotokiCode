clc;  clear; close all;
%% set signed distance function
sdf = @(x) SquareHole(x,2,.5);
%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,5],...
              'NElem',150,...
              'MaxIteration',500,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',.5,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'Type','PlaneStress',...
              'Movie',false,...
              'MovieAxis',[-0.5 2.5 0 5],...
              'MovieCAxis',[-0.2 0.5],...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,3]);

%% material
fem.Material = NeoHookeanMaterial('E',2,'Nu',0.4); 

%% solving
fem.dampedsolve();

function D = SquareHole(P,H,R)
R1 = dRectangle(P,0,H,0,H);
C1 = dCircle(P,H/2,H/2,R);
D = dDiff(R1,C1);
end
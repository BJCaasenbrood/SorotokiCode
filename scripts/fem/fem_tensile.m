clc;  clear; close all;
%% set signed distance function
sdf = @(x) TensileBone(x,8,2,3.5,1,0.75);
%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,10],...
              'NElem',500,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/20,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,8]);

%fem.Material = Ecoflex0030;
%fem.Material = MooneyMaterial('C10',8,'K',160);
fem.Material = NeoHookeanMaterial('E',60,'Nu',0.48); 
%% solving
fem.solve();

function D = TensileBone(P,H,W,T,D,R)
dD = 0.5*(W-D);
dT = 0.5*(H-T);

R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,0,dD,dT,dT+T);
R3 = dRectangle(P,W-dD,W,dT,dT+T);
C1 = dCircle(P,dD-R,dT,R);
C2 = dCircle(P,dD-R,dT+T,R);
C3 = dCircle(P,W-dD+R,dT,R);
C4 = dCircle(P,W-dD+R,dT+T,R);
D0 = dDiff(dDiff(dDiff(R1,R2),C1),C2);
D = dDiff(dDiff(dDiff(D0,R3),C3),C4);
end
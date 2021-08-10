clear; close all; clc;

%% set signed distance function
W = 10;
H = 17.5;
sdf = @(x) Bellow(x,W,H);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',750);
msh = msh.generate();
msh.show(); pause(2);

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,'ResidualNorm',1e-3,'VolumeInfill',0.43,...
              'Penal',4,'VolumetricPressure',true,'FilterRadius',1.25,...
              'Nonlinear',0,'ReflectionPlane',[1,1],'Repeat',[2],...
              'MaxIterationMMA',75,'OptimizationProblem','Compliant','Movie',0);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);

id = fem.FindNodes('Location',[0.01*W,H]);
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,2]);

id = fem.FindNodes('Location',[W,0.01*H]);
fem = fem.AddConstraint('Output',id,[.01,0]);
fem = fem.AddConstraint('Spring',id,[2,0]);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[-1e-2,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0;0,H/2],  6.0);
% fem = fem.initialTopology('Hole',[0,H/3],6.0);
% fem = fem.initialTopology('Hole',[0,H/2],6.0);

%% material
fem.Material = Ecoflex0050;

%% solving
fem.optimize();
fem.show('ISO',0.3);

%% convert topology result to mesh
ISO  = 0.3;
Simp = 0.05;
GrowH = 1.0;
MinH = W/5;
MaxH = H/5;

mshr = fem.exportMesh(ISO,Simp,[GrowH,MinH,MaxH]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/15,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Top'); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindEdges('AllHole');
femr = femr.AddConstraint('Pressure',id,[5*kpa,0]);

id = femr.FindNodes('Bottom');
femr = femr.AddConstraint('Support',id,[1,0]);
femr = femr.AddConstraint('Output',id,[0,0]);

%% assign material to reduced fem
D = 5; % compress. factor (more stable)
femr.Material = Dragonskin30(D);

%% solve final finite-element problem
femr.solve();

function D = Bellow(x,W,H)
R1 = dRectangle(x,0,W,0,H);
C2 = dCircle(x,W,H,W/3);

D = dDiff(R1,C2);
end
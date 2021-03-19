clear; close all; clc;

%% set signed distance function
W = 8;
H = 5;
sdf = @(x) Bellow(x,W,H);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',1200);
msh = msh.generate();
msh.show(); pause(2);

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',0.75,...
              'Nonlinear',false,'TimeStep',1/3,'ReflectionPlane',[1,1],...
              'OptimizationProblem','Compliant','Repeat',[1 2 2],...
              'MaxIterationMMA',50,'Movie',1);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);

id = fem.FindNodes('Location',[0.1*W,H]);
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindNodes('Line',[0.02*W,W,H,H]);
fem = fem.AddConstraint('Spring',id,[0,.1]*1e-1);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[-1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0],1.5);

%% material
fem.Material = Ecoflex0050;

%% solving
fem.optimize();

%% convert topology result to mesh
mshr = fem.exportMesh(0.25,0.07,[1.2,0.25,5]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/15,'FilterRadius',H/15,...
    'MovieAxis',[-75 170 -140 40],'Movie',0,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Top'); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindNodes('Bottom'); 
femr = femr.AddConstraint('Support',id,[1,0]);

id = femr.FindEdges('AllHole');
femr = femr.AddConstraint('Pressure',id,[-0.15*kpa,0]);

%% assign material to reduced fem
D = 25; % compress. factor (more stable)
femr.Material = Ecoflex0050(D);

%% solve final finite-element problem
femr.solve();

function D = Bellow(x,W,H)
R1 = dRectangle(x,0,W,0,H);
C2 = dCircle(x,W,H,H/4);
R2 = dRectangle(x,W/4,W,H*0.9,H);
D = dDiff(dDiff(R1,C2),R2);
end
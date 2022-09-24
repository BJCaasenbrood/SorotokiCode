clear; close all; clc;

P0 = -0.1*kpa;

%% set signed distance function
W = 24;
H = 15;
sdf = @(x) Bellow(x,W,H);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',1200);
msh = msh.generate();
msh.show(); pause(2);

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.225,'Penal',2,'FilterRadius',3,...
              'Nonlinear',1,'TimeStep',1/5,'ReflectionPlane',[1,1],...
              'OptimizationProblem','Compliant','Repeat',[],...%[1,2,2],...
              'MaxIterationMMA',60,'Movie',0,'ChangeMax',0.015);

%% add constraint
fem = fem.addSupport(fem.FindNodes('Bottom'),[0,1]);
fem = fem.addSupport(fem.FindNodes('Left'),[1,0]);

id = fem.FindNodes('Location',[0.1*W,H]);
fem = fem.addOutput(id,[0,sign(P0)]);
fem = fem.addSpring(id,[0,1]);

id = fem.FindNodes('Line',[0.02*W,W,H,H]);
fem = fem.addSpring(id,[0,.1]*1e-1);

id = fem.FindElements('Location',[0,0],1);
fem = fem.addMyocyte(id,[P0,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0],5);

%% material
fem.Material = NeoHookeanMaterial(0.01,0.45);

%% solving
fem.optimize();

%% convert topology result to mesh
mshr = fem.exportMesh(0.25,0.07,[1,2.5,6]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/15,'FilterRadius',H/15,...
    'MovieAxis',[-75 170 -140 40],'Movie',0,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Bottom'); 
femr = femr.addSupport(id,[1,1]);

id = femr.FindNodes('Top'); 
femr = femr.addSupport(id,[1,0]);

id = femr.FindEdges('AllHole');
femr = femr.addPressure(id,P0);

%% assign material to reduced fem
D = 5;
femr.Material = Ecoflex0030(D);

%% solve final finite-element problem
femr.solve();

function D = Bellow(x,W,H)
R1 = dRectangle(x,0,W,0,H);
C2 = dCircle(x,W,H,H/4);
R2 = dRectangle(x,W/4,W,H*0.9,H);
D = dDiff(dDiff(R1,C2),R2);
end
clear; close all; clc;

%% set signed distance function
W = 8;
H = 4;
sdf = @(x) Bellow(x,W,H);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',750);
msh = msh.generate();
msh.show(); pause(2);

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',0.75,...
              'Nonlinear',false,'TimeStep',1/3,'ReflectionPlane',[1,1],...
              'OptimizationProblem','Compliant','Repeat',[1 2],...
              'MaxIterationMMA',65,'Movie',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);

id = fem.FindNodes('Location',[0.01*W,H]);
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,.1]);

id = fem.FindNodes('Line',[0.02*W,W,H,H]);
fem = fem.AddConstraint('Spring',id,[0,.1]*1e-1);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0],1.0);

%% material
fem.Material = Ecoflex0030;

%% solving
fem.optimize();

function D = Bellow(x,W,H)
R1 = dRectangle(x,0,W,0,H);
C2 = dCircle(x,W,H,1.);

D = dDiff(R1,C2);
end
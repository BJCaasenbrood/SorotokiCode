clear; close all; clc;

%% set signed distance function
W = 6;
H = 4;
sdf = @(x) Bellow(x,W,H);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',750);
msh = msh.generate();
msh.show(); pause(2);

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,'ResidualNorm',1e-3,'VolumeInfill',0.3,...
              'Penal',4,'VolumetricPressure',true,'FilterRadius',0.5,...
              'Nonlinear',0,'ReflectionPlane',[1,1],'Repeat',[1 2 2],...
              'MaxIterationMMA',75,'OptimizationProblem','Compliant');

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);

id = fem.FindNodes('Location',[0.01*W,H]);
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[-1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0],1);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();

function D = Bellow(x,W,H)
R1 = dRectangle(x,0,W,0,H);
C2 = dCircle(x,W,H,2);

D = dDiff(R1,C2);
end
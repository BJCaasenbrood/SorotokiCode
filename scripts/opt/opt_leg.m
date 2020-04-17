clear; close all; clc;

%% set signed distance function
sdf = @(x) SoftLeg(x);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,30,0,50],'NElem',1e3);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,...
              'Penal',4,...
              'FilterRadius',5,...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant',...
              'Movie',true);
          
%% set spatial settings
fem = fem.set('Repeat',ones(1,0),'Periodic',[1,0]);

%% add constraint
id = fem.FindNodes('Left');
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[25,0]);
fem = fem.AddConstraint('Support',id,[0,1]);
fem = fem.AddConstraint('Output',id,[-1,0]);
fem = fem.AddConstraint('Spring',id,[1.0,0]);

id = fem.FindElements('Location',[15,30]);
fem = fem.AddConstraint('PressureCell',id,[-1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[15,30],5);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();

function D = SoftLeg(x)
R1 = dRectangle(x,0,30,10,50);
C1 = dCircle(x,15,10,10);
D = dUnion(R1,C1);
end
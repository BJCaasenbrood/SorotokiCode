clear; close all; clc;

%% set signed distance function
sdf = @(x) SoftLeg(x);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,30,0,50],'Quads',650);
msh = msh.generate().show(); pause(3.0)

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',3,'FilterRadius',3,...
              'Nonlinear',false,'MaxIterationMMA',50,...
              'OptimizationProblem','Compliant','ChangeMax',0.01);
          
%% set spatial settings
fem = fem.set('Repeat',ones(5,1),'Periodic',[1/2,0]);

%% add constraint
% id = fem.FindNodes('Location',[30,12]);
% fem = fem.AddConstraint('Support',id,[1,1]);
% id = fem.FindNodes('Location',[0,12]);
% fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Left');
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Right');
fem = fem.AddConstraint('Support',id,[0,1]);

id = fem.FindNodes('Location',[15,0]);
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,1]);
id = fem.FindNodes('Location',[30,25]);
fem = fem.AddConstraint('Output',id,[1,0]);
fem = fem.AddConstraint('Spring',id,[1,0]);

id = fem.FindElements('Location',[10,30]);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[10,30],5);

%% material
fem.Material = Dragonskin10;

%% linear solving
fem.optimize();
fem.show('ISO');

% %% proceed nonlinear
% fem = fem.reset('fem');
% fem = fem.set('Nonlinear',true,...
%               'TimeStep',1/5);
% 
% %% linear solving
% fem.optimize();

function D = SoftLeg(x)
R1 = dRectangle(x,0,30,5,50);
C1 = dRectangle(x,10,20,0,10);
D = dUnion(R1,C1);
end
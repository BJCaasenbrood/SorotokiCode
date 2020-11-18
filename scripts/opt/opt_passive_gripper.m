clr;
%% generate mesh from sdf
sdf = @(x) Gripper(x);

msh = Mesh(sdf,'BdBox',[0,80,0,50],'NElem',1e3);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',7,...
              'Nonlinear',false,'TimeStep',1/3,'ChangeMax',Inf,...
              'OptimizationProblem','Compliance',...
              'MaxIterationMMA',50);

%% set spatial settings
fem = fem.set('ReflectionPlane',[0, 1]);

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

% id = fem.FindNodes('NW'); 
% fem = fem.AddConstraint('Load',id,[1,0]);
% % fem = fem.AddConstraint('Spring',id,[1,0]);

alpha = atan(40/80);
Fn = [-sin(alpha),cos(alpha)];
id = fem.FindNodes('Line',[0 80 0 40]); 
fem = fem.AddConstraint('Load',id,1e-4*Fn);

%% set density
% fem = fem.initialTopology('Hole',[10,30;30,35;50,40],5);

%% material
fem.Material = Ecoflex0050;

%% solving
fem.optimize();
fem.show('ISO');

function Dist = Gripper(P)
R1 = dRectangle(P,0,80,0,50);
R2 = dRectangle(P,70,85,40,55);
C1 = dCircle(P,70,40,10);
L1 = dLine(P,0,80,0,40);

Dist = dIntersect(dUnion(dDiff(R1,R2),C1),L1);
end
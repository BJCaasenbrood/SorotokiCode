clr;
%% generate mesh from sdf
sdf = @(x) Gripper(x);

msh = Mesh(sdf,'BdBox',[0,80,0,50],'NElem',1e3);
msh = msh.generate();

msh.show(); pause(2);

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',7,...
              'Nonlinear',false,'TimeStep',1/3,'ChangeMax',0.05,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',75);

%% set spatial settings
fem = fem.set('ReflectionPlane',[0, 1]);

%% add boundary condition
id = fem.FindNodes('Box',[0 0 40 50]); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Box',[0 0 0 10]); 
fem = fem.AddConstraint('Load',id,[-10,0]);
fem = fem.AddConstraint('Spring',id,[1,0]);

alpha = atan(40/80);
Fn = -[-sin(alpha),cos(alpha)];
id = fem.FindNodes('Line',[0 80 0 40]); 
fem = fem.AddConstraint('Output',id,Fn);
fem = fem.AddConstraint('Spring',id,[0 1]);

%% set density
fem = fem.initialTopology('Hole',[10,30;30,35;50,40],5);

%% material
fem.Material = TPU90;

%% solving
fem.optimize();
fem.show('ISO');

%% convert topology result to mesh
ISO  = 0.3;
Simp = 0.1;
GrowH = 1.1;
MinH = 1.5;
MaxH = 5;

mshr = fem.exportMesh(ISO,Simp,[GrowH,MinH,MaxH]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/12,'FilterRadius',1/15,...
    'Linestyle','none','MovieAxis',[-30 80 0 100],'Movie',1,...
    'PrescribedDisplacement',1);

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Box',[0 0 0 30]); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindNodes('Box',[0 0 70 100]); 
femr = femr.AddConstraint('Support',id,[1,1]);

%fem = femr.AddConstraint('Contact',@(x) SDF(x,20),[0,0]);

id = femr.FindNodes('Box',[0 0 40 60]); 
femr = femr.AddConstraint('Load',id,[-25,0]);

%% assign material to reduced fem
femr.Material = TPU90;

%% solve final finite-element problem
femr.solve();


function Dist = Gripper(P)
R1 = dRectangle(P,0,80,0,50);
R2 = dRectangle(P,70,85,40,55);
C1 = dCircle(P,70,40,10);
L1 = dLine(P,0,80,0,40);

Dist = dIntersect(dUnion(dDiff(R1,R2),C1),L1);
end

function Dist = SDF(P,R)
Dist = dCircle(P,50,50,R);
end
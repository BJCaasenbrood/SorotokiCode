clr;
%% generate mesh from sdf
sdf = @(x) PneuNet(x,20,40,1,20);

msh = Mesh(sdf,'BdBox',[0,20,0,40],'Quads',[25 50]);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',4,...
              'Nonlinear',false,'TimeStep',1/3,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',50,'ChangeMax',0.03,'Movie',0);

%% set spatial settings
fem = fem.set('Periodic',[1/2, 0],'Repeat',ones(7,1));

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);

id = fem.FindElements('Location',[10,25],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[10,25],0.15);

%% material
fem.Material = Ecoflex0030(0.5);

%% solving
fem.optimize();
fem.show('ISO');

function Dist = PneuNet(P,W,H,E,T)
R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,-W/2,E,T,H+H/2);
R3 = dRectangle(P,W-E,W+W/2,T,H+H/2);
C1 = dCircle(P,0,T + 0.5,1);
C2 = dCircle(P,W,T + 0.5,1);
Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end
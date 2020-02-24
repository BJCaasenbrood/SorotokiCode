clr;
%% set signed distance function
sdf = @(x) PneuGrip(x);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-10,13,0,5],'NElem',500);
msh = msh.generate();

%% generate fem model
fem = Fem(msh,'Nonlinear',true,'ReflectionPlane',[0 -1]);
fem = fem.set('VolumeInfill',0.25,'FilterRadius',1.5,...
              'MaxIterationMMA',50,'TimeStep',1/5,...
              'Penal',4,'OptimizationProblem','Compliant');

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,1]);

id = fem.FindNodes('Location',[13,4.5],1);
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindElements('Location',[-10,5],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[-10,5],2);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();

function Dist = PneuGrip(P)
  R1 = dRectangle(P,-10,13,0,5);
  R2 = dRectangle(P,10,15,3,7);
  R3 = dRectangle(P,10,15,0,1);
  C2 = dCircle(P,10,5,2);
  C3 = dCircle(P,10,0,1);
  Dist = dDiff(dDiff(dDiff(dDiff(R1,C2),R2),R3),C3);
end


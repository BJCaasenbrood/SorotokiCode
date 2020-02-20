clr;
%% set signed distance function
sdf = @(x) PneuGrip(x);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-5,13,0,5],'NElem',750);
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh,'Nonlinear',false,'ReflectionPlane',[0 -1]);
fem = fem.set('VolumeInfill',0.3,'FilterRadius',1.5,...
              'Penal',4,'OptimizationProblem','Compliant');

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,1]);

id = fem.FindNodes('Location',[13,4.5],1);
fem = fem.AddConstraint('Output',id,[-1,1]);
fem = fem.AddConstraint('Spring',id,[1,1]);

id = fem.FindElements('Location',[-5,5],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[-5,5],2);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();

function Dist = PneuGrip(P)
  R1 = dRectangle(P,-5,13,0,5);
  R2 = dRectangle(P,10,15,3,7);
  R3 = dRectangle(P,10,15,0,1);
  C2 = dCircle(P,10,5,2);
  C3 = dCircle(P,10,0,1);
  Dist = dDiff(dDiff(dDiff(dDiff(R1,C2),R2),R3),C3);
end


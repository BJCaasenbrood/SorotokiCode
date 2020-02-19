clr;
%% set signed distance function
W = 20;
H = 40;
E = 1;
D = 20;

sdf = @(x) PneuNet(x,W,H,E,D);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,W,0,H],...
              'NElem',1e3,...
              'MaxIteration',50,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.3,...
              'Penal',4,...
              'PrescribedDisplacement',false,...
              'VolumetricPressure',true,...
              'FilterRadius',4,...
              'Periodic',[1/2, 0],...
              'Repeat',[0, 1],...
              'Nonlinear',false,...
              'MaxIterationMMA',150,...
              'Repeat',ones(1,7),...
              'Movie',true,...
              'MovieAxis',[0 160 -45 85],...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindElements('Location',[W/2,H/2],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[W/2,H/2],1);

%% material
fem.Material = Dragonskin10A;

%% solving
figure(101);
fem.optimize();

function Dist = PneuNet(P,W,H,E,D)
R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,-W/2,E,D,H+H/2);
R3 = dRectangle(P,W-E,W+W/2,D,H+H/2);
C1 = dCircle(P,0,D + 0.5,1);
C2 = dCircle(P,W,D + 0.5,1);
Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end


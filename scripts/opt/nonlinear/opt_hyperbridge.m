clr;
%% generate mesh from sdf
sdf = sRectangle(0,80,0,40);

msh = Mesh(sdf,'NElem',3e3);     
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'Nonlinear',0,'Penal',4,'ReflectionPlane',[-1,0],...
              'OptimizationProblem','Compliance','ResidualNorm',Inf,...
              'MaxIterationMMA',30,'ChangeMax',0.1);

fem = fem.set('FilterRadius',1,'VolumeInfill',0.25,'Penal',4);

%% add boundary condition
fem = fem.addSupport(fem.FindNodes('Left'),[1,0]);
fem = fem.addSupport(fem.FindNodes('Right'),[1,1]);
fem = fem.addLoad(fem.FindNodes('NW'),[0,-2e-1]);

%% assign material
fem.Material = NeoHookeanMaterial(1,0.4);

%% set density
%fem = fem.initialTopology('Hole',[1,1;3,1],1);

%% solving
fem.optimize();

%% 
fem.show('ISO')
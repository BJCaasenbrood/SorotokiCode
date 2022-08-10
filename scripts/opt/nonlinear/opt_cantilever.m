clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,120,0,30);

msh = Mesh(sdf,'BdBox',[0,120,0,30],'Quads',[100,20]);     
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'TimeStep',1/30,'Nonlinear',1,...
              'OptimizationProblem','Compliance',...
              'MaxIterationMMA',30,'ChangeMax',0.1);

fem = fem.set('FilterRadius',3,'VolumeInfill',0.3,'Penal',4);

%% add boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('RightMid'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('RightMid'),[0,-8e-3]);

%% assign material
fem.Material = NeoHookeanMaterial(1,0.4);

%% set density
fem = fem.initialTopology('Hole',[1,1;3,1],1);

%% solving
fem.optimize();

%% 
fem.show('ISO')
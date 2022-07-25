clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,2);

msh = Mesh(sdf,'BdBox',[0,5,0,2],'NElem',2000);     
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'TimeStep',1/15,'Nonlinear',1,...
              'OptimizationProblem','Compliance',...
              'MaxIterationMMA',150,'ChangeMax',0.1);

fem = fem.set('FilterRadius',0.1,'VolumeInfill',0.2,...
              'Penal',3,'ReflectionPlane',[-1,0]);

%% add boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',[0,2],1),[0,-1e-3]);

%% assign material
fem.Material = Dragonskin10(15);

%% set density
fem = fem.initialTopology('Hole',[1,1;3,1],1);

%% solving
fem.optimize();

%% 
fem.show('ISO')
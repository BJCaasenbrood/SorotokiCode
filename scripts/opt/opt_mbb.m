clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,2);

msh = Mesh(sdf,'BdBox',[0,5,0,2],'NElem',800);     
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'TimeStep',1/4,'Nonlinear',true,...
              'OptimizationProblem','Compliance');

fem = fem.set('FilterRadius',0.3,'VolumeInfill',0.3,...
              'Penal',4,'ReflectionPlane',[-1,0]);

%% add boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('NW',[4,2],4),[0,-1e-4]);

id = fem.FindElements('Location',[3,0.5],1);
fem = fem.AddConstraint('PressureCell',id,[5e-4,0]);

%% assign material
fem.Material = Dragonskin10A;

%% set density
fem = fem.initialTopology('Hole',[2.5,1],.5);

%% solving
fem.optimize();
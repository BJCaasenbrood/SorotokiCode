clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,2);

msh = Mesh(sdf,'BdBox',[0,5,0,2],'Quads',1000);     
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'TimeStep',1/5,'Nonlinear',true,'FilterRadius',0.3,...
              'OptimizationProblem','Compliance',...
              'MaxIterationMMA',30,'ChangeMax',0.05,'Movie',true);

fem = fem.set('FilterRadius',0.2,'VolumeInfill',0.3,...
              'Penal',4,'ReflectionPlane',[-1,0]);

%% add boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',[4,2],4),[0,-1e-5]);

id = fem.FindElements('Location',[3,0.5],1);

%% assign material
fem.Material = Ecoflex0050;

%% set density
fem = fem.initialTopology('Hole',[2.5,1],.5);

%% solving
fem.optimize();

%% 
fem.show('ISO')
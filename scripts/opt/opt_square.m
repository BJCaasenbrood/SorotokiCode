clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,5,0,10],'NElem',160,'MaxIteration',700,'Triangulate',1);      
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'TimeStep',1/5,'ResidualNorm',1e-3,...
              'FilterRadius',0.1,'VolumeInfill',0.3,...
              'Penal',4,'Nonlinear',1,...
              'ReflectionPlane',[-1 0],...
              'Periodic',[0,1/2],...
              'PrescribedDisplacement',true,...
              'OptimizationProblem','Compliance');

%% add constraint
id = fem.FindNodes('Bottom'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Top');
fem = fem.AddConstraint('Support',id,[1,0]);
fem = fem.AddConstraint('Load',id,[0,3]);

%% material
fem.Material = Dragonskin10A;

%% set density
fem.Density = 0.5*ones(fem.NElem,1);

%% solving
fem.optimize();

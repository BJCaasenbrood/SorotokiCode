%% set signed distance function
sdf = @(x) dRectangle(x,0,8,0,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,8,0,2],'Quads',15);      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'FilterRadius',0.5,...
              'VolumeInfill',0.3,...
              'Penal',4,...
              'OptimizationProblem','Compliance');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[4,2],4); 
fem = fem.AddConstraint('Load',id,[0,-1e-3]);

%% material
fem.Material = Dragonskin10A;

%% set density
%fem = fem.initialTopology('Equidistance',[4,1],0.5);
fem.Density = 0.5*ones(fem.NElem,1);

%% solving
fem.optimize();

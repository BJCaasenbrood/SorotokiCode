%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,2);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,2],...
              'NElem',500,...
              'MaxIteration',150,...
              'ShowMeshing',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'FilterRadius',0.3,...
              'VolumeInfill',0.3,...
              'Penal',4,...
              'OptimizationProblem','Compliance',...
              'ReflectionPlane',[-1,0],...
              'Nonlinear',true);

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,0]);
id = fem.FindNodes('SE'); 
fem = fem.AddConstraint('Support',id,[0,1]);

id = fem.FindNodes('NW',[4,2],5); 
fem = fem.AddConstraint('Load',id,[0,-5e-4]);

%% material
fem.Material = Dragonskin10A;

%% set density
fem = fem.initialTopology('Equidistance',[2,1],.5);

%% solving
fem.optimize();
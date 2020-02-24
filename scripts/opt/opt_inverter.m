%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,5,0,2],'Quads',50);
      
msh = msh.generate();

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/10,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.3,...
              'Penal',1,...
              'PenalMax',4,...
              'PrescribedDisplacement',true,...
              'FilterRadius',0.2,...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Top'); 
fem = fem.AddConstraint('Support',id,[0,1]);

id = fem.FindNodes('Location',[0,0],5);
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[0,2],1);
fem = fem.AddConstraint('Load',id,[1,0]);
fem = fem.AddConstraint('Spring',id,[.1,0]);

id = fem.FindNodes('Location',[5,2],1);
fem = fem.AddConstraint('Output',id,[-1,0]);
fem = fem.AddConstraint('Spring',id,[.1,0]);

%% set density
fem = fem.initialTopology([1,1],.05);

%% material
fem.Material = LinearMaterial('E',10,'Nu',0.49);

%% solving
fem.optimize();



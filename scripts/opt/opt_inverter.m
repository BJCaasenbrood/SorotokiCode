%% set signed distance function
sdf = @(x) dRectangle(x,0,7,0,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,7,0,2],'NElem',1000);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',0.2,...
              'Nonlinear',false,'TimeStep',1/3,'ChangeMax',0.05,...
              'OptimizationProblem','Compliant','MaxIterationMMA',40);


%% set symmetry          
fem = fem.set('ReflectionPlane',[0,-1]);               
%           
%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Top',5),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('NW'),[-1e-3,0]);
fem = fem.AddConstraint('Spring',fem.FindNodes('NW'),[0.1,0]);

fem = fem.AddConstraint('Output',fem.FindNodes('NE'),[1,0]);
fem = fem.AddConstraint('Spring',fem.FindNodes('NE'),[0.1,0]);

%% set density
fem = fem.initialTopology('Hole',[2.5,1],.5);

%% material
fem.Material = Dragonskin30;

%% solving
fem.optimize();



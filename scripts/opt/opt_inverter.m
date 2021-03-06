clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,7,0,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,7,0,2],'Quads',[60 30]);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.4,'Penal',4,'FilterRadius',0.1,...
              'Nonlinear',false,'TimeStep',1/3,'ChangeMax',0.02,...
              'OptimizationProblem','Compliant','MaxIterationMMA',120);


%% set symmetry          
fem = fem.set('ReflectionPlane',[0,-1]);               
%           
%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Top',5),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('NW'),[-1e-3,0]);
fem = fem.AddConstraint('Spring',fem.FindNodes('NW'),[1,0]);

fem = fem.AddConstraint('Output',fem.FindNodes('NE'),[1,0]);
fem = fem.AddConstraint('Spring',fem.FindNodes('NE'),[1,0]);

%% set density
%fem = fem.initialTopology('Hole',[2.5,1],.5);

%% material
fem.Material = TPU90;%LinearMaterial('E',1,'Nu',0.3);

%% solving
fem.optimize();



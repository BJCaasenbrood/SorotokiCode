clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,7,0,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,7,0,2],'Quads',[120 30]);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.4,'Penal',4,'FilterRadius',0.25,...
              'Nonlinear',false,'TimeStep',1/5,'ChangeMax',0.03,...
              'OptimizationProblem','Compliant','MaxIterationMMA',120,...
              'Linestyle','none');

%% set symmetry          
fem = fem.set('ReflectionPlane',[0,-1]);               
%           
%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Top',5),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('NW'),[-0.1,0]);
%fem = fem.AddConstraint('Spring',fem.FindNodes('NW'),[1e-5,0]);

fem = fem.AddConstraint('Output',fem.FindNodes('NE'),[1,0]);
%fem = fem.AddConstraint('Spring',fem.FindNodes('NE'),[1e-5,0]);

%% set density
fem = fem.initialTopology('Hole',[2.5,1],0.4);

%% material
fem.Material = Ecoflex0030(30);

%% solving
fem.optimize();



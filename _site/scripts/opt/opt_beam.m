clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,5,0,2],'NElem',1000);      
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'TimeStep',1/5,'ResidualNorm',1e-3,'FilterRadius',0.3,...
        'VolumeInfill',0.3,'Penal',4,'Nonlinear',0,...
              'OptimizationProblem','Compliance');

%% set symmetry          
fem = fem.set('Periodic',[1, 0],'ReflectionPlane',[-1,0]);          
          
%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('SW'),[0,-1e-5]);

%% material
fem.Material = Ecoflex0050;

%% set density
fem = fem.initialTopology('Hole',[2.5,1],.5);

%% solving
fem.optimize();

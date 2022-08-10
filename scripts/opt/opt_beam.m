clr;
%% generate mesh from sdf
sdf = sRectangle(0,5,0,2);

%% generate mesh
msh = Mesh(sdf,'NElem',10e3);      
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'TimeStep',1/15,'FilterRadius',0.05,...
        'VolumeInfill',0.3,'Penal',4,'Nonlinear',0,...
        'OptimizationProblem','Compliance','MaxIterationMMA',70);

%% set symmetry          
fem = fem.set('Periodic',[1, 0],'ReflectionPlane',[-1,0]);          
          
%% add constraint
fem = fem.addSupport(fem.FindNodes('Left'),[1,0]);
fem = fem.addSupport(fem.FindNodes('Right'),[1,1]);
fem = fem.addLoad(fem.FindNodes('SW'),[0,-1]);

%% material
fem.Material = Elastosil28;

%% set density
fem = fem.initialTopology('Hole',[2.5,1],.5);

%% solving
fem.optimize();

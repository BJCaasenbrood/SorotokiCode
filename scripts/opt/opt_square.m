clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,5,0,10],'Quads',1000);    
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'FilterRadius',0.75,'VolumeInfill',0.2,'MaxIterationMMA',200,...
              'Penal',3,'OptimizationProblem','Compliance','Nonlinear',0);

%% set symmetry          
fem = fem.set('Periodic',[0 1/2],'ReflectionPlane',[-1 0],'Repeat',[2,2]);               
          
%% add constraint
fem = fem.addSupport(fem.FindNodes('Bottom'),[1,1]);
fem = fem.addLoad(fem.FindNodes('Top'),[1e-2,1e-2]);

%% material
fem.Material = NeoHookeanMaterial(1,0.4);

%% set density
fem = fem.initialTopology('Hole',[3,3;3,7;0,2;0,8], 1);

%% solving
fem.optimize();

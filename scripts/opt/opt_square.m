clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,15,0,15);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,15,0,15],'Quads',6000);    
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'FilterRadius',0.3,'VolumeInfill',0.2,'MaxIterationMMA',200,...
              'Penal',3,'OptimizationProblem','Compliance','Nonlinear',0);

%% set symmetry          
fem = fem.set('ReflectionPlane',[0 -1],'Repeat',[2,2]);               
          
%% add constraint
fem = fem.addSupport(fem.FindNodes('Bottom'),[1,1]);
fem = fem.addLoad(fem.FindNodes('Top'),[1e-3,1e-3]);
fem = fem.addGravity();

%% material
fem.Material = NeoHookeanMaterial(1,0.4);

%% set density
%fem = fem.initialTopology('Hole',[3,3;3,7;0,2;0,8], 1);
fem = fem.initialTopology('Random');

%% solving
fem.optimize();

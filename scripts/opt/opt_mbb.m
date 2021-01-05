clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,2);

msh = Mesh(sdf,'BdBox',[0,5,0,2],'NElem',1000);     
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'TimeStep',1/5,'Nonlinear',false,...
              'OptimizationProblem','Compliance',...
              'MaxIterationMMA',100,'ChangeMax',0.03);

fem = fem.set('FilterRadius',0.2,'VolumeInfill',0.4,...
              'Penal',3,'ReflectionPlane',[-1,0]);

%% add boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
%fem = fem.AddConstraint('Load',fem.FindNodes('Location',[0,2],1),[0,-1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Bottom'),[0,-1]);
%% assign material
fem.Material = Dragonskin10(1);

%% set density
fem = fem.initialTopology('Hole',[1,1;3,1],.2);

%% solving
fem.optimize();

%% 
fem.show('ISO')
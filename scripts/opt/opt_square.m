clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,5,0,10],'Quads',1000);    
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'FilterRadius',0.75,'VolumeInfill',0.2,'MaxIterationMMA',100,...
              'Penal',2,'OptimizationProblem','Compliance','Nonlinear',0);

%% set symmetry          
fem = fem.set('Periodic',[0 1/2], 'ReflectionPlane',[-1 0]);               
          
%% add constraint
id = fem.FindNodes('Bottom'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Top');
fem = fem.AddConstraint('Load',id,[3,0]);

%% material
fem.Material = Dragonskin30;

%% set density
fem = fem.initialTopology('Hole',[3,3;3,7;0,2;0,8],1);

%% solving
fem.optimize();

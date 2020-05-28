clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-1,1,-1,1,0,6);
msh = Mesh(sdf,'BdBox',[-1,1,-1,1,0,6],'Hexahedron',[10,10,25]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/3,'PrescribedDisplacement',0,'Nonlinear',0);

fem = fem.set('FilterRadius',0.75,'VolumeInfill',0.3,...
              'Penal',4,'ChangeMax',5);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),so3([0,0,pi/3]));
%fem = fem.AddConstraint('Load',fem.FindNodes('Location',[0,-1,6]),[0,-1e-3,0]);

%% select material
fem.Material =  TPU90();
fem.Density = 0.5*ones(msh.NElem,1);

%% solving
fem.optimize();
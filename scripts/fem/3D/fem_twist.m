clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-1,1,-1,1,0,20);
msh = Mesh(sdf,'BdBox',[-1,1,-1,1,0,20],'Hexahedron',[4,4,30]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/10,'PrescribedDisplacement',false);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,0,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,1e-1,0]);
%fem = fem.AddConstraint('Load',fem.FindNodes('Location',[0,-1,6]),[0,-1e-3,0]);

%% select material
fem.Material =  TPU90();
%fem.Density = 0.5*ones(msh.NElem,1);

%% solving
fem.solve();
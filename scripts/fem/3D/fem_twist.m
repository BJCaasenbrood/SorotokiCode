clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-1,1,-1,1,0,20);
msh = Mesh(sdf,'BdBox',[-1,1,-1,1,0,20],'Hexahedron',[4,4,30]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),so3([0,0,pi/2]));

%% select material
fem.Material =  TPU90();

%% solving
fem.solve();
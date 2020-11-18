clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-3,3,-1,1,0,40);
msh = Mesh(sdf,'BdBox',[-3,3,-1,1,0,40],'Hexahedron',[4,4,30]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/15,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,5,-15]);

%% select material
fem.Material =  TPU90();

%% solving
fem.solve();
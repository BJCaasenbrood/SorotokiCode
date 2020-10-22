clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-3,3,-1,1,0,50);
msh = Mesh(sdf,'BdBox',[-3,3,-1,1,0,50],'Hexahedron',[3,3,15]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'Nonlinear',true,'TimeStep',1/50,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),so3([0,0,1.5*pi]));

%% select material
fem.Material =  TPU90();

%% solving
fem.solve();
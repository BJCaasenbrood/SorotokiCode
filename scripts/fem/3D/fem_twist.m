clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-10,10,-10,10,-10,10);
msh = Mesh(sdf,'BdBox',[-10,10,-10,10,-10,10],'Hexahedron',[10,10,10]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'Nonlinear',true,'TimeStep',1/5,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),so3([0,0,pi]));

%% select material
fem.Material =  Ecoflex0030();

%% solving
fem.solve();
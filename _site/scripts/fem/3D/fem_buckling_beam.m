clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-3,3,-1,1,0,40);
msh = Mesh(sdf,'BdBox',[-3,3,-1,1,0,40],'Hexahedron',[3,3,30]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/30,'PrescribedDisplacement',true,'Movie',true,...
    'MovieAxis',[-6,6,-15,15,0,45]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,0.5,-15]);

%% select material
fem.Material = TPU90();

%% solving
fem.solve();
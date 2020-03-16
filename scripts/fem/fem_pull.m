clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,3,0,1);

msh = Mesh(sdf,'BdBox',[0,3,0,1],'Quads',5,'Triangulate',true);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/25,'PrescribedDisplacement',true);

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[2,0]);

%% assign material
fem.Material = TPU90;

%% solving
fem.solve();
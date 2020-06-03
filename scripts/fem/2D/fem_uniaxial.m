clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,1,0,1);

msh = Mesh(sdf,'BdBox',[0,1,0,1],'Quads',10);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'PrescribedDisplacement',true,...
              'Linestyle','none');

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[4,0]);

%% assign material
fem.Material = TPU90();

%% solve
fem.solve();

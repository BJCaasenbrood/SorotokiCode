clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,1,0,1);

msh = Mesh(sdf,'BdBox',[0,1,0,1],'Quads',3);
msh = msh.generateMesh;

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/25,'PrescribedDisplacement',true);

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[5,0]);

%% assign material
fem.Material = NeoHookeanMaterial('E',3,'Nu',0.4);

%% solving
fem.solve();

fem.show('BC');
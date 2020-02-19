clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,1,0,1);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,1,0,1],'Center',Quads([0,1,0,1],10,2));
msh = msh.generateMesh;

%% generate fem model from mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/25,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[5,0]);

fem.Material = Dragonskin10A;

%% solving
fem.solve();
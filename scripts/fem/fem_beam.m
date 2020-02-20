clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,10,0,1);

msh = Mesh(sdf,'BdBox',[0,10,0,1],'Quads',[25,4]);
msh = msh.generateMesh;

%% add boundary conditions 
fem = Fem(msh,'TimeStep',1/15);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Bottom'),[0,-2e-2]);

%% select material
fem.Material = Dragonskin10A;

%% solving
fem.solve();

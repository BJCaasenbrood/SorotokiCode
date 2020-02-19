clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,1);

msh = Mesh(sdf,'BdBox',[0,5,0,1],'Quads',[25,5]);
msh = msh.generateMesh;

%% add boundary conditions 
fem = Fem(msh,'TimeStep',1/6);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Bottom'),[0,-1e-2]);

%% select material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

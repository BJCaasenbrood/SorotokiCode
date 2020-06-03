clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,10,0,1);

msh = Mesh(sdf,'BdBox',[0,10,0,1],'NElem',150);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/40);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Bottom'),[0,-3e-2]);

%% select material
fem.Material =  Dragonskin10A;

%% solving
fem.solve();

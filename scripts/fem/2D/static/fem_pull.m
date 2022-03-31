clr;
%% generate mesh from sdf
sdf = sRectangle(0,5,0,5);

msh = Mesh(sdf,'BdBox',[0,5,0,5],'Quads',[10 10]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/100);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[0,1]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Right'),[5,0]);

%% select material
fem.Material = Dragonskin10(30);

%% solving
fem.solve();

clr;
%% generate mesh from sdf
sdf = sCube(-3,3,-1,1,0,40);
msh = Mesh(sdf,'Hexahedron',[2,2,50]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/300,'MovieAxis',[-6,6,-25,25,0,45]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0,0]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Top'),[0,1e-1,-25]);

%% select material
fem.Material = Dragonskin10(15); 

%% solving
fem.solve();
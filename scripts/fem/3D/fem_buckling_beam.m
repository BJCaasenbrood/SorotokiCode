clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-3,3,-1,1,0,40);
msh = Mesh(sdf,'BdBox',[-3,3,-1,1,0,40],'Hexahedron',[2,2,30]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'MovieAxis',[-6,6,-15,15,0,45]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0,0]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Top'),[0,1,-15]);

%% select material
fem.Material = Dragonskin10(15); 

%% solving
fem.solve();
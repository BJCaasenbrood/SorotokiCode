clr;
%% generate mesh from sdf
sdf = @(x) dCube(x,-3,3,-3,3,0,20);
msh = Mesh(sdf,'BdBox',[-3,3,-3,3,0,20],'Hexahedron',2*[1,1,2]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'Nonlinear',true,'TimeStep',1/100,...
    'Movie',false,'MovieAxis',[-5 5 -5 5 0 21]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,0,1]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Top'),so3([0,0,0.5*pi]));

%% select material
fem.Material = NeoHookeanMaterial(10,0.33);

%% solving
fem.solve();
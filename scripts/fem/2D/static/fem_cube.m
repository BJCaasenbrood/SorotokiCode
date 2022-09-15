clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,20,0,20);

msh = Mesh(sdf,'BdBox',[0,20,0,20],'Quads',[25 4]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/25);

%% add constraint
fem = fem.addSupport(fem.FindNodes('Left'),[1,1]);
fem = fem.addLoad(fem.FindNodes('Bottom'),[0,-2e-3]);

%% select material
fem.Material = Dragonskin10;

%% solving
fem.solve();

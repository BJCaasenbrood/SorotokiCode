clr;
%% generate mesh from sdf
sdf = sRectangle(0,40,0,2);

msh = Mesh(sdf,'Quads',[50 2]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/120,'SolverPlot',true);

%% add constraint
fem = fem.addSupport(fem.FindNodes('Left'),[1,1]);
fem = fem.addGravity([0,-9.81e3]);

%% select material
fem.Material = NeoHookeanMaterial(0.1,0.4);

%% solving
fem.solve();
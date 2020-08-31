clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,10,0,1);

msh = Mesh(sdf,'BdBox',[0,10,0,1],'Quads',[20 4]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/10);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Bottom'),[0,-5e-3]);

%% select material
fem.Material =  Ecoflex0030;

%% solving
fem.solve();

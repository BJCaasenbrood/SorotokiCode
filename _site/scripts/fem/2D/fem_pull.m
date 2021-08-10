clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,5);

msh = Mesh(sdf,'BdBox',[0,5,0,5],'Quads',[10 10]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/100,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[5,0]);

%% select material
fem.Material = Dragonskin10(100);

%% solving
fem.solve();

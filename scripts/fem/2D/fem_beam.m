clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,40,0,1);

msh = Mesh(sdf,'BdBox',[0,40,0,1],'Quads',[40 4]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/100,'Nonlinear',1);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9810]);

%% select material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

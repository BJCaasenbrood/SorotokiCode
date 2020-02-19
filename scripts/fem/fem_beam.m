clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,1);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,1],'Center',Quads([0,5,0,1],25,5));
msh = msh.generateMesh;

%% add boundary conditions 
fem = Fem(msh);
fem = fem.set('TimeStep',1/50);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Bottom'),[0,-1e-2]);

%% select material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

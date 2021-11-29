clr;
%% settings
P0  = -10*kpa;

%% generate mesh
Simp  = 0.1;
GrowH = 1;
MinH  = 2.5;
MaxH  = 3;

msh = Mesh('Bucklecell.png','BdBox',[0,66.7,0,100],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh.show();

%% re-orient the mesh
%msh = Blender(msh,'Rotate',{'x',90}); msh = msh.show(); pause(2);

%% generate fem model
fem = Fem(msh,'TimeStep',1/175,'Solver','gmres','Linestyle','none');

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);

fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),[P0,0]);

%% assign material
fem.Material = Ecoflex0030(25);

%% solve
fem.solve();

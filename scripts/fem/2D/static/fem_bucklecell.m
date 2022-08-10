clr;
%% settings
P0  = -50*kpa;

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
fem = Fem(msh,'TimeStep',1/30,'Linestyle','none','SolverId',3);

%% add boundary constraint
fem = fem.addSupport(fem.FindNodes('Top'),[1,0]);
fem = fem.addSupport(fem.FindNodes('Bottom'),[1,1]);
fem = fem.addDisplace(fem.FindNodes('Top'),[0,-20]);

%% assign material
fem.Material = NeoHookeanMaterial(1,0.3);

%% solve
fem.solve();

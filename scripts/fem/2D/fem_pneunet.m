clr;
%% generate mesh
Simp  = 0.02;
GrowH = 1;
MinH  = 2;
MaxH  = 3;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh = msh.show();
pause(2);

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/25,'Linestyle','none');

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,15]),[1,1]);

id = fem.FindEdges('AllHole');
fem = fem.AddConstraint('Pressure',id,[25*kpa,0]);

%% assign material
fem.Material = Dragonskin30(25);

%% solve
fem.solve();


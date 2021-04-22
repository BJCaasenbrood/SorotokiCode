clr;
%% generate mesh
Simp  = 0.05;
GrowH = 1;
MinH  = 1;
MaxH  = 1;

msh = Mesh('Octarm_hollow.png','BdBox',[0,120,0,12],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh = msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/10,'Linestyle','none');

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);

id = fem.FindEdges('BoxHole',[2,5]);
fem = fem.AddConstraint('Pressure',id,[0.2*kpa,0]);

%% assign material
fem.Material = Dragonskin10;

%% solve
fem.solve();


clr;

%% generate mesh
Simp = 0.09;
W = 10*4.35;
H = 10*14.35; 

msh = Mesh('SWR_V3.png','BdBox',[0,W,0,H],'SimplifyTol',Simp,...
    'Hmesh',[1.1,1.5,7]);

msh = msh.generate();
msh.show(); 

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/25,'Linestyle','none',...
    'MovieAxis',[-1 45 -1  170],'Movie',1);

%% assign material
fem.Material = TPU90;

%% add boundary constraint
id = fem.FindEdges('Hole',2);
fem = fem.AddConstraint('Support',id{:},[1,1]);

id = fem.FindEdges('Hole',12);
fem = fem.AddConstraint('Support',id{:},[1,0]);

P0 = -25e-3;
% 
id = fem.FindEdges('Hole',2:11);
fem = fem.AddConstraint('Pressure',id,[P0,0]);

%% solve
f = figure(101);
fem.solve();


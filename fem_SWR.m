clr;

%% generate mesh
Simp  = 0.05;
W = 10*4.35;
H = 10*14.35; 
%H = 10*21.39; 

msh = Mesh('SWR_V3.png','BdBox',[0,W,0,H],'SimplifyTol',Simp,...
    'Hmesh',[1.0,1.5,5]);

msh = msh.generate();
msh.show(); 

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/25,'Linestyle','none',...
    'MovieAxis',[-1 45 -1  213.9000],'Movie',0);

%% assign material
fem.Material = TPU90;

%% add boundary constraint
id = fem.FindEdges('Hole',2);
fem = fem.AddConstraint('Support',id{:},[1,1]);

id = fem.FindEdges('Hole',12);
fem = fem.AddConstraint('Support',id{:},[1,0]);
%fem = fem.AddConstraint('Load',id{:},-[0,0.2]);

P0 = -25e-3;
% 
id = fem.FindEdges('Hole',2:11);
fem = fem.AddConstraint('Pressure',id,[P0,0]);

%% solve
f = figure(101);
fem.solve();


clr;
%% experimental meshing
GrowH = 1;
MinH  = 1;
MaxH  = .02;

msh = Mesh('Arm.stl','Hmesh',[GrowH,MinH,MaxH]);
msh = Blender(msh,'Scale',100);
%msh.show();

fem = Fem(msh,'TimeStep',1/25,'Solver','lu');


fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
%fem = fem.AddConstraint('Gravity',[],[0,-9.81e3,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),so3([0,0,pi]));

% assign material
fem.Material = Ecoflex0030(55);

%% solve
f = figure(101);
[~,tmpNode] = fem.solve();
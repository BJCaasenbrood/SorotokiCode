clr;
%% set signed distance function
sdf = sCircle(1);

%% generate mesh
msh = Mesh(sdf,'NElem',150,'Triangulate',1);
msh = msh.generate(); 

%% show mesh
msh.show();
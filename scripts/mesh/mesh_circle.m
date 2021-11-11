clr;
%% set signed distance function
sdf = sCircle(1);

%% generate mesh
msh = Mesh(sdf,'NElem',250);
msh = msh.generate(); 

%% show mesh
msh.show();
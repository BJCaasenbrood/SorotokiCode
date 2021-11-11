clr;
%% set signed distance function
sdf = sRectangle(0,1,0,1);

%% generate mesh
msh = Mesh(sdf,'NElem',250);
msh = msh.generate();
msh.show();

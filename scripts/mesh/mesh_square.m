clr;
%% set signed distance function
sdf = sRectangle(0,1,0,1);

%% generate mesh
msh = Mesh(sdf,'Quads',[20,20]);
msh = msh.generate();
msh.show();

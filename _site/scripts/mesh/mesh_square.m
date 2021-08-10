clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,1,0,1);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-1,1,-1,1],'NElem',250,'Movie',true);
msh = msh.generate();


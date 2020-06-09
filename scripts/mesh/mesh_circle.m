clr;
%% set signed distance function
sdf = @(x) dCircle(x,0,0,1);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-1,1,-1,1],'NElem',150);
msh = msh.generate();
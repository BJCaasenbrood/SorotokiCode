clr;
%% set signed distance function
sdf = @(x) dCircle(x,0,0,1);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-1,1,-1,1],'NElem',250,'ShowProcess',true);
msh = msh.generate(); 

%% show mesh
msh.show();
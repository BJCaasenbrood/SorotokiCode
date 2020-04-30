clr;
%% set signed distance function
R = 1;  % radius

sdf = @(x) dCircle(x,0,0,R);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[-R,R,-R,R],'NElem',150,'ShowMeshing',true);
msh = msh.generate();
clr;
%% set signed distance function
R = 1;  % radius
sdf = @(x) sqrt(x(:,1).^2 + x(:,2).^2) - R^2;

%% generate mesh
msh = Mesh(sdf,'BdBox',[-R,R,-R,R]);
msh.showSDF();
%msh = msh.set('BdBox',[-R,R,-R,R],'NElem',150,'ShowMeshing',true);
%msh = msh.generate();
clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,1,0,1);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,1,0,1],'NElem',150,'ShowMeshing',true);
msh = msh.generate();


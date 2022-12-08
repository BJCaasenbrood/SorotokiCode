clr;
%% set signed distance function
sdf = sRectangle(2).rotate(45) - sCircle(1);
sdf.BdBox = 3*[-1,1,-1,1];

%% generate mesh
msh = Mesh(sdf,'NElem',150,'Movie',true);
msh = msh.generate(); 

%% show mesh
msh.show();
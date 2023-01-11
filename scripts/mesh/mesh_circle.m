clr;
%% set signed distance function
sdf = sCircle(1);
%sdf = sRectangle(2).rotate(45) - sCircle(1);
%sdf.BdBox = 3*[-1,1,-1,1];

%% generate mesh
msh = Mesh(sdf,'NElem',150);
msh = msh.generate(); 

%% show mesh
msh.show();
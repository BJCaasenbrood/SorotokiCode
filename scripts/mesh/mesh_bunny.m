clr;
%% load stl
[f,v] = stlreader('Bunny.stl');
B = round(boxhull(v));

%% set signed distance function
sdf = @(x) dCube(x,B(1),B(2),B(3),B(4),B(5),B(2));

%% generate mesh
msh = Mesh(sdf,'BdBox',B,'Hexahedron',[20 20 60]);
msh = msh.generate();

%% load bunny

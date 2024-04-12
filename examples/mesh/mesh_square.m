clr;
%% make signed distance func.
sdf = sRectangle(10,10);

%% make mesh
msh = Mesh(sdf,'Quads',[9,9]);
msh = msh.generate();

msh.show();
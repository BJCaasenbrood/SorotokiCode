% distance function
file = '+preset/assets/stl/Bunny.stl';

% generate mesh
msh = Mesh(file,'ElementSize',5);
msh = msh.generate;

% show mesh
msh.show();
view(10,10);

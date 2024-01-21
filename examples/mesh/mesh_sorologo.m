% distance function
clr;
file = '+preset/assets/img/sorotoki_bwlogo.png';

% generate mesh
msh = Mesh(file,'BdBox',[0,20,0,20],'ElementSize',.5);
msh = msh.generate();

% show mesh
msh.show();

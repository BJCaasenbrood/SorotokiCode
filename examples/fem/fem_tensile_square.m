clr;
% generate mesh from sdf
sdf = sRectangle(5);
msh = Mesh(sdf,'Quads',250);
msh = msh.generate();

% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/75,'LineStyle','-');

% add boundary conditions
fem = fem.addSupport('sw',[1,1]);
fem = fem.addSupport('bottom',[1,1]);
% fem = fem.addSupport('left',[1,0]);
% fem = fem.addSupport('right',[1,0]);
fem = fem.addDisplace('top',[0,25]);

% assign material
fem = fem.addMaterial( NeoHookean(0.1,0.3) );
% fem = fem.addMaterial( Yeoh );

% solving
fem.solve('MaxIteration',5);
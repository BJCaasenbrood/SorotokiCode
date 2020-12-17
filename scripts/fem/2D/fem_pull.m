clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,3,0,1);

msh = Mesh(sdf,'BdBox',[0,3,0,1],'Quads',[30,8]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/25,'PrescribedDisplacement',true,...
              'Linestyle','none','ColorAxis',[0 1]);

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[10,0]);

%% assign material
fem.Material = Ecoflex0050;

%% solving
fem.solve();
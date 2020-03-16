clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,10],'Quads',10);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Line',[0 4 10 10]),[1,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Line',[0 2 10 10]),[0,-5]);

%% assign material
fem.Material = TPU90;

%% solving
fem.solve();
clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,10],'Quads',100);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/75,'Nonlinear',true,'Linestyle','none',...
              'PrescribedDisplacement',true,'ColorAxis',[0 0.04]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Line',[0 10 10 10]),[0,-4]);

%% assign material
fem.Material = Ecoflex0050(5);

%% solving
fem.solve();
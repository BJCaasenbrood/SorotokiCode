clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,10],'Quads',300);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/75,'Nonlinear',true,'Linestyle','-',...
          'ColorAxis',[0 0.4],'SigmoidFactor',-0.5);

%% add constraint
fem = fem.addSupport(fem.FindNodes('Top'),[1,0]);
fem = fem.addSupport(fem.FindNodes('Bottom'),[0,1]);
fem = fem.addSupport(fem.FindNodes('Left'),[1,0]);
fem = fem.addDisplace(fem.FindNodes('Box',[0 5 10 10]),[0,-5]);

%% assign material
fem.Material = Ecoflex0050(12);

%% solving
fem.solve();
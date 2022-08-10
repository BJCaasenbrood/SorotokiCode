clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,20,0,1);

msh = Mesh(sdf,'BdBox',[0,20,0,1],'NElem',120);
msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/100);

%% add constraint
fem = fem.addSupport(fem.FindNodes('Left'),[1,1]);
fem = fem.addSupport(fem.FindNodes('Right'),[0,1]);
fem = fem.addDisplace(fem.FindNodes('Right'),[-2,0]);

%% assign material
fem.Material = Dragonskin30();

%% solving
fem.solve();
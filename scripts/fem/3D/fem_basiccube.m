clr;
%% meshing
msh = Mesh('BaseCube.stl','Hmesh',[1,1,2]);

%% fem
fem = Fem(msh,'TimeStep',1/30,'Linestyle','-');
fem.Material = Ecoflex0030(1);

%%
fem = fem.addSupport(fem.FindNodes('Location',[-1,-1,-1]),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Bottom'),[0,0,1]);
fem = fem.addDisplace(fem.FindNodes('Top'),[0,0,1]);

%%
fem.solve();
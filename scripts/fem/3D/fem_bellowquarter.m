clr;
%% meshing
msh = Mesh('BellowQuarter.stl','Hmesh',[1,1,3]);
msh.generate();

%% fem
fem = Fem(msh);
fem.Material = NeoHookeanMaterial();

%%
fem = fem.addSupport(fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.addDisplace(fem.FindNodes('Top'),[0,0,12]);

%%
fem.solve();
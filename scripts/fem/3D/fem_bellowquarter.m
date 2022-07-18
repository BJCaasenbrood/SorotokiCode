clr;
%% meshing
msh = Mesh('BellowQuarter.stl','Hmesh',[1,2,3]);
msh.generate();

%% fem
fem = Fem(msh);
fem.Material = NeoHookeanMaterial();

%%
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Top'),[0,0,12]);

%%
fem.solve();
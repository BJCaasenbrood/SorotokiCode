clr;
%% meshing
msh = Mesh('TensileBone.stl','Hmesh',[1,3,6]);

%% fem
fem = Fem(msh,'TimeStep',1/120,'Linestyle','-');
fem.Material = Ecoflex0030(10);

%%
fem = fem.AddConstraint('Support',fem.FindNodes('Location',[-9.5,-2,0]),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,1,0]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Top'),[0,0,-40]);

%%
fem.solve();
clr;
%% meshing
sdf = sCube(1);
obj = Gmodel(sdf,'Quality',40);
obj.export('tmp.stl');

%%
msh = Mesh('tmp.stl','Hmesh',[1,0.5,0.75]);
delete('tmp.stl');

%% fem
fem = Fem(msh,'TimeStep',1/30,'Linestyle','-',...
    'MovieAxis',[-2,2,-2,2,-2,2]);

fem.Material = Ecoflex0030(25);

%%
fem = fem.AddConstraint('Support',fem.FindNodes('Location',[-1,-1,-1]),[1,1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,0,1]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Top'),[0,0,-1]);

%%
fem.solve();
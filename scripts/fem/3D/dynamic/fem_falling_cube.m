clr;
%% meshing
msh = Mesh('BaseCube.stl','Hmesh',[1,1,2]);

%% fem
fem = Fem(msh,'TimeStep',1/300,'TimeEnd',5,'Linestyle','-');
fem.Material = Ecoflex0030(1);

%%
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,0]);
fem = fem.AddConstraint('Gravity',[],[0,0,-9.81e3]);
fem = fem.AddConstraint('Contact',@(x) SDF(x),[0,0,0]);

%%
fem.simulate();


function d = SDF(x)
sdf = sCube(-4,4,-4,4,-2,-1);
d = sdf.eval(x);
d = d(:,end);
end
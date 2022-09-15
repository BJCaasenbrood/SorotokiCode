clr;
%% meshing
msh = Mesh('TensileBone.stl','Hmesh',[1,5,6]);
msh = msh.generate();

%% fem
fem = Fem(msh,'TimeStep',1/350,'Linestyle','-');
fem.Material = Dragonskin10(1);

%%
%fem = fem.addSupport(fem.FindNodes('Location',[-9.5,-2,0]),[1,1,1]);
fem = fem.addSpring(fem.FindNodes('Left'),[1,1,1]);
%fem = fem.addSupport(fem.FindNodes('Left'),[0,1,0]);
%fem = fem.addSupport(fem.FindNodes('Top'),[1,1,0]);
fem = fem.addLoad(fem.FindNodes('Right'),[0,0,1]);
%fem = fem.addDisplace(fem.FindNodes('Top'),so3([0,0,2*pi]));
%fem = fem.addGravity([0,0,1]);
%%
fem.simulate();
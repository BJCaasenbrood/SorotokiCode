clr;
%% meshing
msh = Mesh('TensileBone.stl','Hmesh',[1,3,6]);
msh = msh.generate();
msh = Blender(msh,'Rotate',{'x',90});
msh = Blender(msh,'Rotate',{'y',90});

%% fem
fem = Fem(msh,'TimeStep',1/50,'Linestyle','-');
fem.Material = NeoHookeanMaterial(1,0.3);

%%
L = fem.FindNodes('Front');
R = fem.FindNodes('Back');

fem = fem.addSupport(L,[1,1,1]);
fem = fem.addDisplace(R,[0,100,0]);

%%
fem.solve();
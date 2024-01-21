clr;

msh = preset.mesh.pneunet();

fem = Fem(msh,'TimeStep',1/60,'BdBox',[-50 120 -100 20]);
fem = fem.addMaterial(NeoHookean(0.2, 0.33));
fem = fem.addMaterial(NeoHookean(5., 0.33));

bottomlayer = msh.findElements('box',[0,120,0,1.75]);
fem = fem.setMaterial(bottomlayer,2);

fem = fem.addPressure('allhole', 10 * 1e-3);
fem = fem.addSupport('left',[1,1]);

fem = fem.addContact(sCircle(20,[25,-30]));

fem.options.LineStyle = '-';
fem.options.isNonlinear = 1;
fem.options.Display = @plt;

fem = fem.solve;

function plt(Fem)
    cla;
    showMaterialsFem(Fem);
    colorbar off;
    showContactFem(Fem);
end

clr;
msh = Mesh(sRectangle(10, 2),'Quads',[30,10]);
msh = msh.generate();

fem = Fem(msh,'SpatialFilterRadius',0.3);

fem = fem.addMaterial(NeoHookean);
fem = fem.addSupport('se',[1, 1]);
fem = fem.addSupport('sw',[1, 1]);
fem = fem.addLoad('bottommid',[0,1]);

fem.options.isNonlinear = false;
fem.options.LineStyle   = 'none';
fem.options.Display     = @plt;

fem = fem.optimize();

function plt(Fem)
    cla;
    figure(101);
    showInfillFem(Fem);
end



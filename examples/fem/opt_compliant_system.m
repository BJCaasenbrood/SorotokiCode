clr;
msh = Mesh(sRectangle(10, 10),'Quads',[55,55]);
msh = msh.generate();

fem = Fem(msh,'SpatialFilterRadius',0.3);

fem = fem.addMaterial(NeoHookean);
fem = fem.addSupport('nw',[1, 1]);
fem = fem.addSupport('sw',[1, 1]);
fem = fem.addLoad('leftmid',[-1,0]);
fem = fem.addOutput('rightmid',[1,0]);

fem.topology.Type = 'compliant';
fem.options.isNonlinear = false;
fem.options.LineStyle = 'none';
fem.options.Display = @plt;
fem.topology.MaxChange = 0.15;
fem = fem.optimize;

function plt(Fem)
    cla;
    showInfillFem(Fem);
end


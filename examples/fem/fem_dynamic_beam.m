clr;
sdf = sRectangle(0, 100, 0, 5);
msh = Mesh(sdf,'Quads',[50,3]);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/1e3);
fem = fem.addMaterial(NeoHookean(.1,0.3));
fem = fem.addSupport('left', [1, 1]);
fem = fem.addContact(sCircle(7,[40,-25]));
fem = fem.addGravity();

fem.options.LineStyle = 'none';
% fem.options.ColorMap  = cmap_redblue;
fem.options.Display = @plt;

fem = fem.simulate('MaxIteration',30);

function plt(Fem)
    clf;
    showVonMisesFem(Fem);
    showContactFem(Fem);
end
clr;
sdf = sRectangle(0, 100, 0, 5);
msh = Mesh(sdf,'Quads',[50,3]);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/500);
fem = fem.addMaterial(NeoHookean(.01,0.45));
fem = fem.addSupport('left', [1, 1]);
fem = fem.addContact(sCircle(7,[40,-25]));
fem = fem.addGravity();

fem.options.Display = @plt;
fem.BdBox = [-2, 120, -80, 10];

fem = fem.simulate('MaxIteration',4);

function plt(Fem)
    clf;
    Fem.show;
    showContactFem(Fem);
end
clr;
sdf = sRectangle(0, 100, 0, 5);
msh = Mesh(sdf,'Quads',[50,3]);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/500);
fem = fem.addMaterial(NeoHookean(.05,0.45));
fem = fem.addSupport('left', [1, 1]);
fem = fem.addContact(sCircle(7,[40,-25]));
fem = fem.addGravity();

fem.options.Display = @plt;
fem.BdBox = [-2, 120, -80, 10];

fem = fem.simulate;

function plt(Fem)
    clf;
    Fem.show;
    if isfield(Fem.system,'Contact')
        showContactFem(Fem);
    end
end
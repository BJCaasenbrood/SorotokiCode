clr;
% Create a mesh object
msh = preset.mesh.suzumori_walker('n',1250);

% Create a Fem object
fem = Fem(msh,'TimeStep',1/3333,'LineStyle','none');

% Add a Neo-Hookean material
fem.addMaterial(NeoHookean(1.2, 0.49));

% Set the density and viscous damping
fem.materials.Material{1}.params.Rho  = 3600e-12;
fem.materials.Material{1}.params.Zeta = 1.75;

% Set the contact parameters
fem.materials.Material{1}.contact.NormalReaction  = 0.5;
fem.materials.Material{1}.contact.TangentFriction = 0.5;

% Add gravity
gvec = -9800*[sind(-30); cosd(-30)];
fem = fem.addGravity(gvec);

flr = sLine(60,0,-3,-3);
fem = fem.addContact(flr);

% Set the display function
fem.options.Display = @(x) plt(x,p);
fem.solver.MaxIteration = 10;

% solve dynamic simulation
fem = fem.simulate();

function plt(fem,p)
    clf;
    fem.showVonMises;
    fem.showContact;
    view(30,90)
end



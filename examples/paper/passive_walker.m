clr;
% Create a mesh object
msh = preset.mesh.suzumori_walker('n',1250);

% Create a Fem object
fem = Fem(msh,'TimeStep',1/3333,'LineStyle','none');

% Add a Neo-Hookean material
fem.addMaterial(NeoHookean(1.2, 0.45));

% Set the density and viscous damping
fem.materials.Material{1}.params.Rho  = 3600e-12;
fem.materials.Material{1}.params.Zeta = 1.75;

% Set the contact parameters
fem.materials.Material{1}.contact.TangentReaction = 1.0;
fem.materials.Material{1}.contact.NormalReaction  = 0.1;
fem.materials.Material{1}.contact.ContactFriction = 1.0;

% Add gravity
gvec = -9800*[sind(-30); cosd(-30)];
fem = fem.addGravity(gvec);

flr = sLine(60,0,-3,-3);
fem = fem.addContact(flr);

% Set the display function
fem.options.Display = @(x) plt(x);
fem.solver.MaxIteration = 20;

% solve dynamic simulation
fem = fem.simulate();

function plt(fem)
    clf;
    fem.showVonMises;
    fem.showContact;
    % view(30,90)
end



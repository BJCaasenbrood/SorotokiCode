clr;

% mesh
sdf = sCircle(5,[18,25]) - sCircle(3.5,[18,25]);
msh = Mesh(sdf,'NElem',55);
msh = msh.generate();

% Create a Fem object
fem = Fem(msh,'TimeStep',1/750);
% Add a Neo-Hookean material
fem = fem.addMaterial(NeoHookean(0.01,0.3));

% Set the contact parameters
fem.materials.Material{1}.contact.NormalReaction  = 1;
fem.materials.Material{1}.contact.TangentReaction = 1;
fem.materials.Material{1}.contact.ContactFriction = 0.3;

% Add gravity
fem = fem.addGravity();
fem = fem.addContact(SDF(15),[0,0]);

% Set the display function
fem.options.Display = @plt;

% simulate
fem = fem.simulate();

% Signed distance environment
function D = SDF(W)
    S1 = sLine(4,-W,-10,-5);
    S2 = sLine(W,4,0,-10);
    S3 = sLine(-W,-W,-1,1);
    D = S1 + S2 + S3;
end

% assigns the plot function to the fem object
function plt(fem)
    cla;
    fem.showVonMises;
    fem.showContact;
    ylim([-15,30]);
end

function fem = suzumori_walker(varargin)
% Create a parser for user inputs
p = inputParser;
% Add optional inputs and default values
addOptional(p,'theta',-30);
addOptional(p,'n',1500);
addOptional(p,'d',120);
addOptional(p,'dt',1/322);

addOptional(p,'gravity',1);
addOptional(p,'contact',1);
addOptional(p,'hanging',0);

% Parse the inputs
parse(p,varargin{:});

% Create a mesh object
msh = preset.mesh.suzumori_walker('n',p.Results.n);

% Create a Fem object
fem = Fem(msh,'TimeStep',p.Results.dt,'LineStyle','none');

% Add a Neo-Hookean material
fem.addMaterial(NeoHookean(1.2, 0.49));

% Set the density and viscous damping
fem.materials.Material{1}.params.Rho  = 3600e-12;
fem.materials.Material{1}.params.Zeta = 1.75;

% Set the contact parameters
fem.materials.Material{1}.contact.NormalReaction  = 0.5;
fem.materials.Material{1}.contact.TangentFriction = 0.5;

% Add gravity
if p.Results.gravity
    gvec = -9800*[sind(p.Results.theta); cosd(p.Results.theta)];
    fem = fem.addGravity(gvec);
end

% Add a contact plane
if p.Results.contact
    flr = sLine(p.Results.d,0,-3,-3);
    fem = fem.addContact(flr);
end

% Set the initial configuration
if p.Results.hanging
    fem.addSupport(fem.findNodes('Location',[40,13],20),[1,1]);
end

% Set the display function
fem.options.Display = @(x) plt(x,p);
fem.solver.MaxIteration = 10;

% log = Log('Stamp',[]);

% log.info(['Material = ', fem.materials.Material{1}.Type]);
% log.info(['Timestep = ' num2str(1/fem.solver.TimeStep,3) ' hz']);
% log.info(['Gravity  = ', num2str(9800,3) ' mm/s^3']);

% assigns the plot function to the fem object
function plt(fem,p)
    clf;
    fem.showVonMises;
    
    if p.Results.contact
        fem.showContact;
    end

    view(4.5,90)
end

end
function fem = multigait_crawler(varargin)
% Create a parser for user inputs
p = inputParser;
% Add optional inputs and default values
addOptional(p,'n',0.75);
addOptional(p,'d',250);
addOptional(p,'dt',1/750);
addOptional(p,'control',0);

% Parse the inputs
parse(p,varargin{:});

% Create a mesh object
msh = preset.mesh.multigait_crawler('n',p.Results.n);

% Create a Fem object
fem = Fem(msh,'TimeStep',p.Results.dt);
% Add a Neo-Hookean material
fem.addMaterial(NeoHookean(0.2, 0.4));
fem.addMaterial(NeoHookean(1.0, 0.4));

% set bottom layer to material 2
BotLayer = fem.findElements('Box',[0,150,0,2]);
fem = fem.setMaterial(BotLayer,2);

% Add gravity
fem = fem.addGravity();

% Add a contact plane
flr = sLine(p.Results.d,-10,-.01,-.01);
fem = fem.addContact(flr);

% Set the display function
fem.options.Display = @plt;
% fem.solver.MaxIteration = 6;

if p.Results.control
    fem = fem.addPressure(fem.findEdges('BoxHole',[0 50 0 15]),...
        @(t) pset(t,1));
    fem = fem.addPressure(fem.findEdges('BoxHole',[50 100 0 15]),...
        @(t) 0.75*pset(t,2));
    fem = fem.addPressure(fem.findEdges('BoxHole',[100 150 0 15]),...
        @(t) pset(t,3));
end

fem.BdBox = [-5,250,-5, 50];

end

function y = pset(t,id)
w  = 2*pi;
P0 = 15 * 1e-3;
t  = mod(t,.75);
y  = P0*tsin(w*t - (id-1)*pi/5).*sigmoid(3*w*t);
end

% assigns the plot function to the fem object
function plt(fem)
cla;
fem.showVonMises;
fem.showContact;
axis equal;
end
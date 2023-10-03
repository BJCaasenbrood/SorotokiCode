function fem = bouncing_ball(varargin)
% Create a parser for user inputs

p = inputParser;
% Add optional inputs and default values
addOptional(p,'n',1500);
addOptional(p,'d',250);
addOptional(p,'dt',1/1250);

% Parse the inputs
parse(p,varargin{:});

% mesh
sdf = sCircle(5,[18,25]) - sCircle(3.5,[18,25]);
msh = Mesh(sdf,'NElem',55);
msh = msh.generate();

% Create a Fem object
fem = Fem(msh,'TimeStep',p.Results.dt);
% Add a Neo-Hookean material
fem = fem.addMaterial(NeoHookean(0.01,0.3));

% Set the contact parameters
% fem.materials.Material{1}.contact.NormalReaction  = 1.0;
% fem.materials.Material{1}.contact.TangentFriction = 0.5;

% Add gravity
fem = fem.addGravity();
fem = fem.addContact(SDF(15),[0,0]);

% Set the display function
fem.options.Display = @plt;

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
    end
end
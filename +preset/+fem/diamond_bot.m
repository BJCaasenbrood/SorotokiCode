function fem = diamond_bot(varargin)
    % Create a parser for user inputs
    p = inputParser;
    % Add optional inputs and default values
    addOptional(p,'dt',1/25);
    addOptional(p,'n',5);
    addOptional(p,'gravity',0);
    addOptional(p,'contact',0);

    % Parse the inputs
    parse(p,varargin{:});
    
    % Create a mesh object
    msh = preset.mesh.diamond_bot('n',p.Results.n);
    
    % Create a Fem object
    fem = Fem(msh,'TimeStep',p.Results.dt);

    % Add a Neo-Hookean material
    fem = fem.addMaterial(NeoHookean(1., 0.33));

    % add support on left side
    fem = fem.addSupport('bottom',[1,1]);
        
    % Add gravity
    if p.Results.gravity
        gvec = -9800*[sind(p.Results.theta); cosd(p.Results.theta)];
        fem = fem.addGravity(gvec);
    end

    fem.options.LineStyle = 'none';
    fem.options.Display = @(x) plt(x,p);
    % fem.solver.MaxIteration = 10;

    % assigns the plot function to the fem object
    function plt(fem, p)
        cla;
        fem.showVonMises;
        
        if p.Results.contact
            fem.showContact;
        end
    end
    
    end
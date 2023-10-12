function fem = pneunet(varargin)
    % Create a parser for user inputs
    p = inputParser;
    % Add optional inputs and default values
    addOptional(p,'dt',1/75);
    addOptional(p,'gravity',0);
    addOptional(p,'contact',0);

    % Parse the inputs
    parse(p,varargin{:});
    
    % Create a mesh object
    msh = preset.mesh.pneunet;
    
    % Create a Fem object
    fem = Fem(msh,'TimeStep',p.Results.dt);

    % Add a Neo-Hookean material
    fem = fem.addMaterial(NeoHookean(1,  0.33));
    % fem = fem.addMaterial(NeoHookean(5., 0.33));

    % % add inextensible bottom layer
    % bottomlayer = msh.findElements('box',[0,120,0,1.75]);
    % fem = fem.setMaterial(bottomlayer,2);

    % add support on left side
    fem = fem.addSupport('left',[1,1]);

    % Add gravity
    if p.Results.gravity
        gvec = -9800*[sind(p.Results.theta); cosd(p.Results.theta)];
        fem = fem.addGravity(gvec);
    end
    
    if p.Results.contact
       fem = fem.addContact(sCircle(20, [25,-30]));
    end

    % Add pressure loads to all holes
    fem = fem.addPressure('allhole', @(t) 10 * 1e-3 * smoothstep(t) );

    fem.options.LineStyle = 'none';
    fem.options.Display = @(x) plt(x,p);

    % assigns the plot function to the fem object
    function plt(fem, p)
        cla;
        fem.show;
        
        if p.Results.contact
            fem.showContact;
        end
        axis equal;
    end
end
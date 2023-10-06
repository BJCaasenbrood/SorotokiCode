function shp = octarm(varargin)

    p = inputParser;
    addOptional(p,'n',30);
    addOptional(p,'m',10);
    parse(p,varargin{:});

    % try
        mat = NeoHookean(1,0.49);
        shp = Shapes(chebyspace(p.Results.n, p.Results.m),[0,p.Results.m,0,0,0,0],...
            'Material',mat,'Length',165,'Texture',matcap_softmath);
        
        shp = shp.setRadius([8,8,0.8]);
        % shp.system.Drag        = 100e-12;
        shp.solver.TimeStep    = 1/60;
        shp.solver.TimeHorizon = 10;
        
        shp = shp.rebuild();        
    % end
end
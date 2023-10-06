function shp = threebellow(varargin)

    p = inputParser;
    addOptional(p,'n',30);
    addOptional(p,'m',1);
    addOptional(p,'Texture',matcap_diffuse(0.75));
    addOptional(p,'isHQ',false);
    parse(p,varargin{:});

    Shading = 'Vertex';
    currentDir = fileparts(mfilename('fullpath'));
    if p.Results.isHQ
        
        filePath =fullfile(currentDir, '..', 'assets', 'stl', 'threebellow_hq.stl');   
        Shading = 'Face';
    else
        filePath =fullfile(currentDir, '..', 'assets', 'stl', 'threebellow.stl');   
    end

    % try
        obj = Gmodel(filePath,'Texture',p.Results.Texture,'Shading',Shading);
        mat = NeoHookean(1,0.49);
        shp = Shapes(chebyspace(p.Results.n, p.Results.m),[0,1,1,0,0,0],...
            'Material',mat,'Length',60);
        
        shp = addGmodel(shp,obj);
        shp = shp.setRadius([12,12,0]);
        shp.solver.TimeStep    = 1/60;
        shp.solver.TimeHorizon = 10;
        
        shp = shp.rebuild();        
    % end
end
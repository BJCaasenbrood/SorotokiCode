function shp = threebellow_finger(varargin)

    p = inputParser;
    addOptional(p,'n',10);
    addOptional(p,'m',1);
    addOptional(p,'phi',0);
    addOptional(p,'q0',1e-2);
    addOptional(p,'isHQ',false);
    addOptional(p,'Texture',matcap_diffuse(0.75));
    
    parse(p,varargin{:});
    
    Shading = 'Vertex';
    currentDir = fileparts(mfilename('fullpath'));
    if p.Results.isHQ
        filePath =fullfile(currentDir, '..', 'assets', 'stl', 'threebellow_finger_hqq.stl');    
    else
        filePath =fullfile(currentDir, '..', 'assets', 'stl', 'threebellow_finger.stl');   
    end

    % try
        obj = Gmodel(filePath,'Texture',p.Results.Texture,'Shading',Shading);
        mat = NeoHookean(1,0.49);
        shp = Shapes(chebyspace(p.Results.n, p.Results.m),[0,0,1,0,0,0],...
            'Material',mat,'Length',45,'Texture',matcap_softmath);
        
        deg = (30) * (pi/180);
        shp = shp.setBase(SE3(rotx(p.Results.phi * 2*pi/3)) * ...
            SE3(rotx(pi/2) * roty(deg) * rotx(-pi/2),[16.5,16.5,0]));

        shp = addGmodel(shp,obj);
        shp = shp.setRadius([12,12,0]);
        shp.solver.TimeStep = 1/60;
        shp.solver.TimeHorizon = 10;
        
        shp = shp.rebuild();        

        % shp.solver.sol.x(1) = p.Results.q0;
    % end
end
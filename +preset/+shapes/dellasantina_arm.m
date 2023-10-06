function shp = dellasantina_arm(varargin)

    p = inputParser;
    addOptional(p,'n',50);
    addOptional(p,'dt',1/120);
    parse(p,varargin{:});

    % try
        M  = 6;
        Y  = pccspace(30,M);
        dL = 63;

        shp = Shapes(Y,[0,0,6,0,0,0], 'Length', M * dL,...
            'VolumetricContact', false, ...
            'Texture', matcap_diffuse(0.45));

        shp.Material = NeoHookean(5,0.33);
        shp.Material.params.Rho = 50e-12;
        shp.Material.params.Zeta = 0.1;
        
        shp.Material.contact.ContactFriction = 0;
        shp.Material.contact.NormalReaction  = 0.1;
        shp.Material.contact.TangentReaction = 0;

        shp = shp.setRadius([22, 15, 0.75]);
        shp = shp.setBase(SE3(rotz(pi/2),[0;0;0]));

        % t   = linspace(0,1,100);
        % shp = shp.setRamp(smoothfall(10*t-9.35));
        shp = shp.setRamp(1e-3);

        [sdf, V] = environmentSdf;
        shp = shp.addContact(sdf);
        shp = shp.rebuild();

        shp.solver.TimeHorizon = Inf;
        shp.solver.TimeStep = p.Results.dt;
        shp.solver.MaxIteration = 10;

        figure(101);
        set_colororder;
        fplot(V,'Color',col(1),'LineW',6);
        hold on;
    %end
end

function [sdf, V] = environmentSdf
    currentDir = fileparts(mfilename('fullpath'));
    functionPath = fullfile(currentDir,'assets');
    addpath(functionPath);

    V = wavy_wall();
    sdf = sPolyline(flipud(V));
    sdf.BdBox = [-245 0 200 390];
end
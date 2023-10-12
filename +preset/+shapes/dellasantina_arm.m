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

    V = wavy_wall;
    sdf = sPolyline(flipud(V));
    sdf.BdBox = [-245 0 200 390];
end

function V = wavy_wall

    V = [-265.2150  243.2586
    -260.1783  244.5286
    -255.1586  245.8644
    -250.1518  247.2478
    -245.1565  248.6721
    -240.1728  250.1363
    -235.1984  251.6318
    -230.2343  253.1615
    -225.2898  254.7529
    -220.3453  256.3446
    -215.3798  257.8695
    -210.4890  259.6143
    -206.0999  262.3194
    -203.0810  266.5246
    -200.9680  271.2658
    -199.2287  276.1601
    -197.4676  281.0464
    -195.2980  285.7617
    -192.2283  289.9308
    -188.0058  292.9216
    -183.2098  294.9029
    -178.2285  296.3725
    -173.2108  297.7157
    -168.2578  299.2762
    -163.5573  301.4684
    -159.5735  304.7699
    -156.6419  309.0433
    -154.5403  313.7895
    -152.8211  318.6909
    -151.0210  323.5624
    -148.6736  328.1891
    -145.3821  332.1872
    -141.0401  335.0032
    -136.1869  336.8441
    -131.1942  338.2761
    -126.1889  339.6645
    -121.2770  341.3473
    -116.7051  343.7889
    -113.2476  347.6172
    -110.8075  352.1992
    -108.5849  356.8940
    -106.4485  361.6285
    -104.4381  366.4180
    -102.5043  371.2390
    -100.6350  376.0853
     -98.8249  380.9541
     -97.0652  385.8413
     -95.3524  390.7452
     -93.6856  395.6649
     -92.0713  400.6020];
end
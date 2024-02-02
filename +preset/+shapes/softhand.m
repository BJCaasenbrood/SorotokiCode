function [SHAPES, GRIPPER] = softhand(varargin)

clr;
p = inputParser;
addOptional(p,'n',25);
addOptional(p,'dt',1/60);
addOptional(p,'contact',1);
parse(p,varargin{:});

currentDir = fileparts(mfilename('fullpath'));
stlPath =fullfile(currentDir, '..', 'assets', 'stl', 'softhand_mount.stl');    

    obj = Gmodel(stlPath,'Shading','Face','Texture',matcap_egg);
    figure(101);
    view(105,20);
    obj.render();

    L = [90;90;90;90;65];
    
    GRIPPER = obj;
    clearvars obj
    
    G = initialBase;
    Y = chebyspace(p.Results.n,2);
    SHAPES = cell(4,1);

    for ii = 1:5
        shp = Shapes(Y,[0,2,0,0,0,0],'Length',L(ii),...
          'Texture', matcap_diffuse(0.45) );
        
        shp = shp.setInputMap( @(x) [1; 0] );
        shp.Material = NeoHookean(1.5, 0.3);
        shp.Material.params.Rho = 150e-12;
        shp.Material.params.Zeta = 0.05;
        shp = shp.setRadius([10,10,0,.75]);
        shp = shp.addGravity();
        shp = shp.setBase(G{ii});

        phi = @(k) -(k-1)*pi/6;
        shp = shp.rebuild();

         shp = shp.addControl( @(x) Control(x, phi(ii)) );

        shp.solver.TimeStep = 1/30;
        shp.solver.MaxIteration = 2;
        shp.solver.TimeHorizon = Inf;

        SHAPES{ii} = showRenderShapes(shp);
    end
    axis([-20 20 -10 130 -100 100]);

end

% cell function containing SE(3) bases for soft fingers
function G = initialBase
G = cell(5);
for ii = 1:4
    G{ii} = SE3(eye(3),[0,22*(ii-1),0]) * SE3(roty(-pi/2));
end
    G{5} = SE3(rotx(-pi/4)*rotz(0),[6,83,-38]) * SE3(roty(-pi/2));
end

% controller function for periodic swinging
function u  = Control(shp, phi)
    t = shp.solver.Time;
    omega = pi;
    Pset  = 350;
    u = sign(sin(omega*t - phi(1)));   
    u = smoothstep(t)*clamp(u*Pset,5,Pset);
end

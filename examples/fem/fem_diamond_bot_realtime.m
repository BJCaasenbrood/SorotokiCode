clr;
% generate fem model from diamond bot
msh = preset.mesh.diamond_bot('n',1e3);
fem = Fem(msh,'isNonlinear',false);
fem = fem.addSupport('bottom',[1,1]);

% build input and output matrices
[B, G, C, Ia] = buildInputOutput(fem);
fem.system.InputMap = @(x) B;

% assigning materials
mat = NeoHookean(2.5, 0.35);
mat.params.Zeta = 2;
fem = fem.addMaterial(mat);

fem = fem.compute(); 
fem.BdBox = [-40 200 0 160];

% control  loop
[U, w, ki, dt] = deal([0;0], 2*pi, 7.5 * mat.getModulus, 1/120);
[uout, pout, pdout] = deal([]);

while fem.solver.Time < Inf

    t = fem.solver.Time;

    % compute error;
    e = C.' * fem.solver.sol.x(Ia) - Pd();
    H = C.' * (fem.system.Tangent \ G);

    % I-action controller
    U = U - dt * ki * (H \ e);

    fem.system.fControl = U;
    fem = fem.update(dt);

    plt(fem, Pd());
end

%% input and output function
function [B,G,C,Ia] = buildInputOutput(fem)
    B = zeros(fem.Mesh.NNode*fem.Dim,2);
    C = zeros(fem.Mesh.NNode*fem.Dim,2);
    qa1 = fem.findNodes('location',[5,50]) * fem.Dim - 1;
    qa2 = fem.findNodes('location',[165,50]) * fem.Dim - 1;
    px  = fem.findNodes('location',[85,95]) * fem.Dim - 1;
    py  = fem.findNodes('location',[85,95]) * fem.Dim;
    
    B(qa1,1) = 1;
    B(qa2,2) = -1;
    C(px,1)  = 1;
    C(py,2)  = 1;
    
    Ia = fem.system.Ia;
    G = B(Ia,:);
    C = C(Ia,:);  
end

%% plotting function
function plt(fem, Pd)

    cla;
    xmin = 30; xmax = 140;
    ymin = 70; ymax = 140;

    patch('Vertices',[xmin,ymin; xmax,ymin; xmax,ymax; xmin,ymax; xmin,ymin],...
        'Faces',[1,2,3,4,5],'FaceColor',[1,1,1] * 0.875,'LineStyle','none'); 

    fem.show();
    fplot(Pd' + [85,95],'k.','MarkerSize',45,'Color',col(2),'LineW',2);

    if fem.solver.Time < 1/120
        axis equal;
    end
    
    axis( [-40 200 0 160] );
    drawnow;
end

function out = Pd
    xmin = 30; xmax = 140;
    ymin = 70; ymax = 140;
    x = get(gca, 'CurrentPoint');
    if (x(1,1) < -abs(2*xmin))
        out = [0; 10];
    else
        out = [0;0];
        out(1) = clamp(x(1,1)-85,xmin-85, xmax-85);
        out(2) = clamp(x(1,2)-95,ymin-95, ymax-95);
        %out = [x(1,1); x(1,2)] - [85; 95];
    end
end

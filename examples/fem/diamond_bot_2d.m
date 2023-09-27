clr;
% generate fem model from diamond bot
msh = preset.mesh.diamond_bot('n',250);
fem = Fem(msh);
fem = fem.addSupport('bottom',[1,1]);

% build input and output matrices
[B, G, C, Ia] = buildInputOutput(fem);
fem.system.InputMap = @(x) B;

% assigning materials
mat = NeoHookean(2.0, 0.45);
mat.params.Zeta = 0.01;
fem = fem.addMaterial(mat);

fem = fem.compute();    % precompute
fem = fem.set('MaxIteration',10);

% control  loop
[U, w, kt, dt] = deal([0;0], 2*pi, 10 * mat.getModulus, 1/120);
[uout, pout, pdout] = deal([]);

while fem.solver.Time < 2

    t = fem.solver.Time;
    Pd = [25 * cos(2 * w * t); 25 * sin(w * t)];

    % compute error;
    p = C.' * fem.solver.sol.x(Ia);
    e = p - Pd;
    H = C.' * (fem.system.Tangent \ G);

    % I-action controller
    U = U - dt * ki * (H \ e);
    uout  = [uout; U.'];
    pout  = [pout; p.'];
    pdout = [pdout; Pd.'];

    fem.system.fControl = U;
    fem = fem.update(dt);

    plt(fem, Pd);
end

figure(102);
plot(uout(:,1),'LineW',2); hold on;
plot(uout(:,2),'LineW',2); 

figure(103);
subplot(2,1,1);
plot(pout(:,1),'LineW',2); hold on;
plot(pdout(:,1),'LineW',2); 
subplot(2,1,2);
plot(pout(:,2),'LineW',2); hold on;
plot(pdout(:,2),'LineW',2); 

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
    fem.showVonMises;
    clim([-1e7,1e7]);
    Nds = fem.Mesh.Node + meshfield(fem, fem.solver.sol.x);
    id = fem.findNodes('location',[85,95]);
    fplot(Nds(id,:),'.','MarkerSize',60,'Color',col(1));
    fplot(Pd' + [85,95],'k.','MarkerSize',30,'Color',col(2));
    
    if fem.solver.Time < 1/120
         axis([-40 200 0 160]);
    end

    drawnow;
end

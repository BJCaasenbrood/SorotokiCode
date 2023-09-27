clr;
% setting up the simulation
Xd1 = [-140.6; 353.7; 0];
Xd2 = [-230.7; 257.6; 0];

% load preset shapes
shp = preset.shapes.dellasantina_arm;
fig = figure(101);

madeContact = false;

% simulation loop
while shp.solver.Time < shp.solver.TimeHorizon
    
    % contact enabler
    if norm(shp.system.fContact) < eps && ~madeContact
        Xd = [-140.6;353.7;0];
    else % if contact, then interp to new position
        madeContact = true;
        Xd1 = [-140.6;353.7;0];
        Xd2 = [-257.6;240.7;0];
        Xd = lerp(Xd1, Xd2, shp.solver.Time/3.5);
    end

    shp = shp.addControl(@(shp) Control(shp, Xd));

    shp = shp.update();
    shp = showRenderShapes(shp);
    
    background();
    axis([-275 25 0 450]);
    drawnow;
    
    % set fig.Name to simulation time located at fem.solver.Time
    fig.Name = sprintf('Time: %1.2f', shp.solver.Time);
end

function tau = Control(shp, Xd)

    t = shp.solver.Time;
    M  = shp.system.Mass;
    C  = shp.system.Coriolis;
    fe = shp.system.fElastic;
    J  = shp.system.Jacobian(4:6,:,end);
    dJ = shp.system.dJacobian(4:6,:,end);
    x  = shp.system.Backbone(1:3,4,end);

    q  = shp.solver.sol.x;
    dq = shp.solver.sol.dx;
    nq = numel(q);
    
    Kc = 1e-3 * smoothstep(2.5*t);
    Dc = 0.1  * Kc;
    
    % compute lambda
    Mi  = inv(M);
    JT  = J.';
    A   = J*Mi*JT + 1e-4 * eye(3);
    lam = A\(eye(3));
    
    % compute JB+
    W = Mi*JT;
    Jbp = W*lam;
    
    % compute eta
    eta = lam*(J*Mi*C - 0*dJ);
    
    % impedance controller
    tau = J.'*Jbp.'*(fe) + J.'*eta*(eye(nq) ...
         - Jbp*J)*dq + J.'*(Kc*(Xd - x) - Dc*J*dq);
end
clr;
%% settings
enableContact = false; % enables sphere contact
Xd = [50,0,20]; % setpoint coordinates
M  = 6; % number of basis modes 

%% Example: PD-Jacobian feedback control for soft manipulator
shp = Shapes(chebyspace(20,M),[0,M,0,0,0,0],...
    'Material',NeoHookean(0.01,0.3));

shp = shp.addGravity([0;0;-500]);
shp = shp.setBase('+x');
shp = shp.setRadius([5,5,0.75]);

shp = shp.rebuild();
shp.solver.TimeStep = 1/30;
shp.solver.TimeHorizon = 15;
shp.solver.MaxIteration = 5;

obj = sSphere(2,Xd);
obj = Gmodel(obj,'Texture',matcap_diffuse(0.82));
obj.render();

if enableContact
    sdf = sSphere(12.5,[50,0,-25]);
    shp = shp.addContact(sdf);
    con = Gmodel(sdf,'Texture',matcap_diffuse(0.2));
    con.render();
end

shp.system.Controller = @(x) Control(x,Xd);
shp = shp.simulate();

T = shp.solver.sol.tout;
Y = shp.solver.sol.yout;

%%
function tau = Control(shp,Xd)
    J = shp.system.Jacobian(4:6,:,end);
    ve = shp.system.Velocity(4:6,:,end);
    fe = shp.system.fElastic;
    fg = shp.system.fBody;
    x  = shp.system.Backbone(1:3,4,end);

    [kd, kp] = deal(1e-5);
    
    tau = fe + fg + kp * J.'*(Xd(:) - x) - kd * J.' * ve;
    % axis tight;
    % view(30,30);
end
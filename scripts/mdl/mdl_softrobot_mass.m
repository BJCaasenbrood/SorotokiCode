clr;

%%
J = eye(3)*0.01;
m = eye(3)*0.02;
L = 120;

%%
rgb = RigidBody('Mtt',blkdiag(J,m));
rgb = rgb.setInitialSE3(SE3(eye(3),[0,0,L]));
rgb = rgb.setGravity();

rgb.Gmodel = Gmodel(sCylinder(0,0,-10,10,10).');

%%
Y   = chebyspace(50,1);
shp = Shapes(Y,[0,1,0,0,0,0],'Length',L);
%shp = shp.setBase(roty(pi/2));
shp = shp.rebuild();

mdl = Model(shp,'TimeEnd',3,'TimeStep',1/250);
mdl = mdl.addSystem(rgb);

mdl.Controller = @(x) fnc(x);

mdl = mdl.simulate();
%%
figure(101);
t = mdl.Log.t;
view(30,30);
axis([0 120,-20,20,-60,60]*2);
rgb.Gmodel.bake.render();

for ii = 1:fps(t,120):numel(t)
    X = mdl.Log.x(ii,:);
    
    q = X(1);

    rgb = rgb.render(X(3:end));
    shp = shp.render(q);
    
    view(30,30);
    axis([0 120,-20,20,-60,60]*2);
end


%%
function u = fnc(mdl)
t = mdl.t;
u = zeros(mdl.NIn,1);
u(1) = 500*sin(4*pi*t);

m   = mdl.Systems{2}.getMass();
ge  = mdl.Systems{1}.Log.FK.g(:,:,end);
eta = mdl.Systems{1}.Log.FK.eta(:,:,end);
gc  = mdl.Systems{2}.Log.g;
etac = mdl.Systems{2}.Log.eta;
R   = gc(1:3,1:3);

eta(1:3)  = eta(1:3);
etac(1:3) = etac(1:3);

Kp = kron(diag([10,20]),eye(3));
Xi = logmapSE3(gc\ge);
Fu = Kp*tmapSE3(Xi)*wedge(Xi) - 2*Kp*admapinv(gc)*(etac - eta)*mdl.TimeStep;

if t < 1.245
    u(2:end) = Fu;
    u(5:end) = u(5:end) + R.'*(m*[0,0,9810].');
end

end
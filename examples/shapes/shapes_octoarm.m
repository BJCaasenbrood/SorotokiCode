clr;
M = 12;

mat = NeoHookean(0.01,0.4);
mat.params.Zeta = 0.1;
shp = Shapes(pccspace(50,M), [0,M,0,0,0,0], 'Length',175);

shp = shp.setMaterial(mat);
shp = shp.setRadius([8,8,0.85]);

shp = shp.addDrag(1200e-12, 0.49, 0.01);
shp.options.isRampCompensation = true;
shp.options.Quality = 60;

shp.solver.TimeStep = 1/30;
shp.solver.TimeHorizon = 50;
shp.solver.MaxIteration = 2;

shp.system.Controller = @(x) Control(x);

figure(101); 
view(30,30);
axis([-120 120 -8 8 -5 150]);
drawnow;

jj = 1;
while shp.solver.Time < shp.solver.TimeHorizon
    shp = shp.update();
    shp = showRenderShapes(shp);
end 

function tau = Control(shp)
    t = shp.solver.Time;
    tau = zeros(shp.NJoint,1);
    
    tau(1:5) = -0.75 * sin(3*t);
end
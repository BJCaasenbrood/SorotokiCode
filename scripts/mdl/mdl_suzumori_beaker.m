clr;
%%
sdf = sCylinder(0,0,0,90,38) + sCylinder(0,0,90,94,41);
sdf = sdf.rotate('y',-0.05);
sdf = sdf.translate([0,0,-120]);

hld = Gmodel('Suzumori_Holder_Grippers.stl','Shading','Face');
hld.Texture = (metalclean)*0.65;
rod = Gmodel('Suzumori_Beam.stl','Shading','Face');
rod.Texture = (metalclean)*0.9;
bkr = Gmodel('Beaker.stl');
bkr = Blender(bkr,'Rotate',{'z',-30});
bkr = Blender(bkr,'Rotate',{'y',-2.85});
bkr = Blender(bkr,'Translate',{'z',-120});

%%
N = 20;
M = 3;

Y = chebyspace(N,M);
shp = Shapes(Y,[0,0,M,0,0,0],...
    'Length',80,...
    'Texture',egg*0.85,...
    'xia0',[0,0,0,1,0,0]);

shp = shp.setRadius(9);
shp.g0 = SE3(roty(pi/2),[0;37.5;5]);
%shp.g0 = SE3(roty(pi),[0;37.5;5]);

shp.Material = NeoHookeanMaterial(0.2,0.3);
shp.Material.Rho = 2000e-12;
shp.Material.Zeta = .75;
shp.Material.Cfr  = 1;
shp.Material.Rr   = 1;
shp.Contact       = sdf;

t   = linspace(0,1,100);
shp = shp.setRamp(smoothfall(10*t-9.35));
shp = shp.rebuild();
shp.X0            = [0.03;0;0;0;0;0];

for ii = 1:4
    theta = (2*pi/4)*(ii-1);
    SHP{ii} = shp;
    SHP{ii}.g0 = SE3(rotz(theta))*SHP{ii}.g0;
end

%%
mdl = Model(SHP{1},'TimeEnd',2,'TimeStep',1/60);

for ii = 2:4
    mdl = mdl.addSystem(SHP{ii});
end

mdl.Controller = @(x) Controller(x);

%%
mdl = mdl.simulate;

%%
figure(102); cla;
Q = mdl.Log.x;
for jj = 1
    Q = mdl.Log.x;
    X = Q(:,(1:3) + (jj-1)*6);
    dX = Q(:,(4:6) + (jj-1)*6);
end

subplot(1,2,1);
plot(mdl.Log.t,X);

subplot(1,2,2);
plot(mdl.Log.t,dX);

%%
rod.render();
hld.render(); 
bkr.render();

axis tight;
view(-10,-2)
drawnow;
rod.update();
hld.update();
bkr.update();
%%
t = mdl.Log.t;

for ii = 1:fps(t,30):numel(t)
      
   for jj = 1:4
     Q = mdl.Log.x(ii,:);
     XX = Q((1:3) + (jj-1)*6);
     mdl.Systems{jj} = mdl.Systems{jj}.render(XX,1);
   end
   
   timetitle(t(ii));
   background('w');
   axis([-50 50  -50  50 -120  80]);

end

%%
function u = Controller(mdl)
    A = 120;
    t = mdl.t;
    u = zeros(12,1);    
    u(1)  = A*(smoothfall(t) - 0.5);
    u(4)  = A*(smoothfall(t) - 0.5);
    u(7)  = A*(smoothfall(t) - 0.5);
    u(10) = A*(smoothfall(t) - 0.5);
end



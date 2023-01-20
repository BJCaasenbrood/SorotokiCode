clr;
%% parameters
M  = 0.012;        % Mass screw (kg)
R  = 5;            % Radius screw (mm)
I  = 0.5*M*R^2;    % Inertia screw
mu = 2;            % Friction screw
E0 = 1;            % Young's modulus soft robot (MPa)
Nu0 = 0.3;         % Poisson ratio (-)
Rho = 2000e-12;    % Density (kg/mm2)
Zeta = 0.3;        % Damping coeff.
Rr = 0.5;          % Reaction force coeff. contact
Fr = 0.5;          % Normal friction coeff. contact

%%
sdf = sCylinder(0,0,0,44,28).translate([0,0,-80]);

hld = Gmodel('Suzumori_Holder_Grippers.stl',...
    'Shading','Face','Texture',(metalclean)*0.65);
rod = Gmodel('Suzumori_Beam.stl',...
    'Shading','Face','Texture',(metalclean)*0.9);

bkr = Gmodel('M10.stl','Shading','Face');
bkr = Blender(bkr,'Translate',{'z',-70});
bkr = Blender(bkr,'Fix');

num = Gmodel('M10_numbers_reverse.stl',...
    'Shading','Face','Texture',0.2*diffuse(0.5));
num = Blender(num,'Translate',{'z',-70});
num = Blender(num,'Fix');

%%
N = 30;
M = 3;

Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,M,0,0,0],'xia0',[0,0,0,1,0,0],...
    'Length',80,'Texture',egg*0.85);

shp    = shp.setRadius(9);
shp.g0 = SE3(roty(pi/2),[0;37.5;5]);

shp.Material = NeoHookeanMaterial(E0,Nu0);
shp.Material.Rho  = Rho;
shp.Material.Zeta = Zeta;
shp.Material.Cfr  = Fr;
shp.Material.Rr   = Rr;
shp.Contact       = sdf;

shp = shp.rebuild();

for ii = 1:4
    theta = (2*pi/4)*(ii-1);
    SHP{ii} = shp;
    SHP{ii}.g0 = SE3(rotz(theta))*SHP{ii}.g0;
end

%%
mdl = Model(SHP{1},'TimeEnd',15,'TimeStep',1/60,'MaxIteration',10);

for ii = 2:4
    mdl = mdl.addSystem(SHP{ii});
end

%% screw model
f = @(x,u,t) [x(2); ...
              (-mu*sign(x(2)) + u)/I];

sys = StateSpace(f,[0,0],'NInput',1);
mdl = mdl.addSystem(sys);

%% controller
mdl.Controller = @(x) Controller(x);
mdl = mdl.simulate;

%%
figure(102); cla;
Q = mdl.Log.x;
for jj = 1
    Q  = mdl.Log.x;
    X  = Q(:,(1:3) + (jj-1)*6);
    dX = Q(:,(4:6) + (jj-1)*6);
end

subplot(1,2,1);
plot(mdl.Log.t,X);

subplot(1,2,2);
plot(mdl.Log.t,dX);

figure(103);
plot(mdl.Log.t,wrapToPi((Q(:,end-1))));

%%
%sdf.show;
rod.render();
hld.render(); 
bkr.render();
num.render();

axis tight;
view(0,5);
drawnow;
rod.update();
hld.update();
bkr.update();
%%
t = mdl.Log.t;
view(10,5);
Pitch = 0.00;

for ii = 1:fps(t,9):numel(t)
      
   for jj = 1:4
     Q = mdl.Log.x(ii,:);
     XX = Q((1:6) + (jj-1)*12);
     mdl.Systems{jj} = mdl.Systems{jj}.render(XX,1);
   end
   
   bkr = bkr.reset; num = num.reset;
   Theta = rad2deg(Q(end-1));

   bkr = Blender(bkr,'Rotate',{'z',-Theta});
   num = Blender(num,'Rotate',{'z',-Theta});
   bkr.update; 
   num.update;

   timetitle(t(ii));
   background('w');
   view(10,5);
   axis([-50 50  -50  50 -160  60]);

end

%%
function u = Controller(mdl)
    W = 2*pi;
    A = 15*8;
    B = 20*8;

    t = mdl.t;
    u = zeros(sum(vertcat(mdl.NInput{:})),1);    

    y1 = -clamp(sin(-W*t),0,1);
    y2 = -clamp(sin(-W*t-pi/2),-.25,1);

    u(4:6:end) = A*y1;
    u(1:6:end) = B*y2;

    u = u*smoothstep(2*t);
    
    tau = 0;
    for ii = 1:4
        tauW = zeros(3,1);
        sys = mdl.Systems{ii};
        V  = sys.Log.Con.Node;
        Ft = sys.Log.Con.Wt;
        
        for jj = 1:size(V,1)
            r = V(jj,:) - [0,0,-70];
            tauW = tauW + hat(r)*Ft(jj,:).';
        end

        tau = tau + tauW(3);
    end

    u(end) = tau;

end



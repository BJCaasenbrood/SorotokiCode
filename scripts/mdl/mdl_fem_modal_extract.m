clr;
%% pressure settings
P0  = 35*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;

MinH  = 1.5;
MaxH  = 3;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],...
    'SimplifyTol',Simp,'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/120,...
              'BdBox',[0,120,-80,20],'TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,1,0,10]),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time));

%% add output nodes
id  = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = NeoHookeanMaterial(1,0.25); 
fem.Material.Zeta = 0.4;

%% solve
fem = fem.simulate(); 

%% shape function reconstruction
shp = Shapes(fem,[0,5,0,2,0,0],'NNode',120,...
     'Sdf',CrossSection,'Center',[7.5;0;0]);

%shp = Shapes(fem,[0,5,0,2,0,0],'NNode',200);
shp = shp.reference([0,2],[120,2]);
shp = shp.rebuild();
shp = shp.reconstruct();

Y = shp.get('PODR');
Z = shp.get('PODQ');

figure(101);
%clf;
for ii = 1:4
   subplot(1,2,1);
   plot(shp.Sigma,Y(:,ii),'LineW',3,'Color',col(ii)); 
   hold on;
   
   subplot(1,2,2);
   plot(shp.Sigma,Z(:,ii),'LineW',3,'Color',col(ii)); 
   hold on;
end

%% model
shp.E  = 5;
shp.Nu = 0.1;
shp.Zeta = 0.5;

shp = shp.rebuild();

mdl = Model(shp,'Tstep',1/60,'Tsim',15);

mdl.tau = @(M) Controller(M);

% simulate system
mdl = mdl.simulate(); 

FPS = 120;

figure(103);
h = [];
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    delete(h);
    
    p = shp.FK(mdl.Log.q(ii,:));
    h = plot(p(:,1),p(:,3),'LineW',3);
    axis equal;
    %view(30,30);
    axis([-40 120 -120 40]);
    drawnow;
    
end

%% setup controller
function tau = Controller(mdl)
n = numel(mdl.Log.q);
t = mdl.Log.t;

tau = zeros(n,1);
tau(1) = 350*smoothstep(t)*sin(2*t) - 100;
tau(2) = 0.1*(350*smoothstep(t)*sin(2*t) - 100);

end

function sdf = CrossSection
H = 20/2;
W = 15;

R1 = sRectangle(0,W,0,H);
R2 = sRectangle(2,W-2,5,H-3);
sdf = R1 - R2;
end

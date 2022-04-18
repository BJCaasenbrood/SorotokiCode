clr;
%% pressure settings
P0  = 20*kpa;

%% generate mesh
MinH  = 3;
MaxH  = 4;

msh = Mesh('Pneunet.png','BdBox',[0,107,0,16],...
    'SimplifyTol',0.02,'Hmesh',[1,MinH,MaxH]);

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
shp = Shapes(fem,[0,2,0,1,0,0],'NNode',120,...
     'Sdf',CrossSection,'Center',[7.5;0;0]);

shp = shp.reference([0,2],[107,2]);
shp = shp.rebuild();
shp = shp.reconstruct();

X = shp.get('PODR');
Y = shp.get('PODQ');

figure(101); clf;

for ii = 1:2
   subplot(1,2,1);
   plot(shp.Sigma,X(:,ii),'LineW',3,'Color',col(ii)); 
   hold on;
   
   subplot(1,2,2);
   plot(shp.Sigma,Y(:,ii),'LineW',3,'Color',col(ii)); 
   hold on;
end

Z = shp.get('POD');
gm = ones(shp.NNode,1);
G = trapz(shp.Sigma.',gm.*Z).';

%% model
shp.E    = 5;
shp.Nu   = 0.1;
shp.Zeta = 0.05;

shp = shp.rebuild();

mdl = Model(shp,'Tstep',1/60,'Tsim',8);

mdl.tau = @(M) Controller(M,G);

mdl = mdl.simulate(); 

FPS = 120;

figure(103);
h = [];
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    delete(h);
    
    p = shp.FK(mdl.Log.q(ii,:));
    h = plot(p(:,1),p(:,3),'LineW',3);
    axis equal;

    axis([-40 120 -120 40]);
    drawnow;
end

%% setup controller
function tau = Controller(mdl,G)
n = numel(mdl.Log.q);
t = mdl.Log.t;

P0 = 1500;
tau = -G*P0*sigmoid(max(t-3,0));
end

function sdf = CrossSection
H = 20/2;
W = 15;

R1 = sRectangle(0,W,0,H);
R2 = sRectangle(2,W-2,5,H-3);
sdf = R1 - R2;
end

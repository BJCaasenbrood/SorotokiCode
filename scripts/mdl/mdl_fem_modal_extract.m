clr;
%% pressure settings
P0  = 50*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;

MinH  = 0.8586*2;
MaxH  = 1.5*2;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/120,'BdBox',[0,120,-80,20],'TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,20]),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,5,0,0]),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,10,20,20]),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[10.9,10.9,9,20]),[1,1]);
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
shp = Shapes(fem,[0,5,0,3,0,0],'NNode',300,...
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
shp = shp.rebuild();

mdl = Model(shp,'Tstep',1/60,'Tsim',5);

%% simulate system
mdl = mdl.simulate(); 

%% 
figure(100);
subplot(2,1,1);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);

subplot(2,1,2);
plot(mdl.Log.t,mdl.Log.q(:,M+1:2*M),'LineW',2);
colororder(col);



function sdf = CrossSection
H = 20/2;
W = 15;

R1 = sRectangle(0,W,0,H);
R2 = sRectangle(2,W-2,5,H-3);
sdf = R1;
end

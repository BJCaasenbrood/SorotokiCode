clr;
%% simulation settings
P = 15*kpa;
W = 120;
H = 20;

%% finite element settings
Simp  = 0.02;
GrowH = 1;
MinH  = 2;
MaxH  = 3;

%% generate mesh
msh = Mesh('Pneunet.png','BdBox',[0,W,0,H],...
    'SimplifyTol',Simp,'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

figure(101);
subplot(2,1,1); imshow('Pneunet.png');
subplot(2,1,2); msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/120,'TimeEnd',2,'SelfCollision',1);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('Hole'),P);
%fem = fem.AddConstraint('Gravity',[],[0,-9810]);

%% assign material
fem.Material = Ecoflex0030(300);
%fem.Material = NeoHookeanMaterial(0.5,0.4);

%% solve
close all;
fem.solve();

%% movie
t = fem.Log.t; 
close all;
figure(105);

for ii = 1:fps(t,60):numel(t)
    N = fem.Log.Node{ii};

    fem.set('Node',N);
    fem.show('Field',fem.Log.Stress{ii});
    axis([-60 130 -120 30]);
    background(gitpage);
    
    drawnow();
end
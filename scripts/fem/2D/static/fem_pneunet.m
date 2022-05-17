
clr;
%% pressure settings
P0  = 30*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;
MinH  = 3;
MaxH  = 4;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

figure(101);
subplot(2,1,1); imshow('Pneunet.png');
subplot(2,1,2); msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/50,'BdBox',[0,120,-80,20],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',2,'SolverPlot',1);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),P0);

%% assign material
fem.Material = Dragonskin30(25);

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
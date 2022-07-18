clr;
%% simulation settings
P = 15*kpa;
W = 120;
H = 20;

%% finite element settings
Simp  = 0.02;
GrowH = 1;
MinH  = 1.5;
MaxH  = 2.5;

%% generate mesh
msh = Mesh('Pneunet.png','BdBox',[0,W,0,H],...
    'SimplifyTol',Simp,'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

figure(101);
subplot(2,1,1); imshow('Pneunet.png');
subplot(2,1,2); msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/120,'TimeEnd',1,'Linestyle','none');

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('Hole'),P);

%% assign material
fem.Material = Dragonskin10();

%% solve
fem.solve();

%% movie
figure(105);
t = fem.Log.t; 

Pos = [];

for ii = 1:fps(t,60):numel(t)
    Pos(:,:,ii) = fem.Log.Node{ii}(fem.FindNodes('Bottom'),:);
    
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    axis([-30 130 -100 30]);   
    background();
    drawnow();
end

figure(106);
C = barney(-numel(t));
for ii = 1:fps(t,20):numel(t)
    plot(Pos(:,1,ii),Pos(:,2,ii),...
        'LineW',2,'Color',C(ii,:)); hold on;
end
axis equal; axis([-30 130 -100 30]);   


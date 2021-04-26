clr;
%% generate mesh
Simp  = 0.02;
GrowH = 1;
MinH  = 2;
MaxH  = 3;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh = msh.show();
pause(2);

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/15,'Linestyle','none');

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,15]),[1,1]);

id = fem.FindEdges('AllHole');
fem = fem.AddConstraint('Pressure',id,[40*kpa,0]);

%% add output nodes
id = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Output',fem.FindNodes('Bottom'),[0,0]);

%% assign material
fem.Material = Dragonskin30();

%% solve
fem.solve();

%% post-processing
Ux = fem.Log{2,2};
Uy = fem.Log{3,2};
N0 = fem.get('Node0');

figure(103); cla;
for ii = 1:1:size(Ux,2)
    Nx = N0(id,1) + Ux(:,ii);
    Ny = N0(id,2) + Uy(:,ii);
    plot(Nx,Ny,'Linewidth',2,...
        'Color',col(4,ii/(1.05*size(Ux,2))));
    hold on;
end

axis([-20 120 -90 10]); axis equal;
xaxis('$x$-dimension','mm');
yaxis('$y$-dimension','mm');
set(gca,'linewidth',1.5)
grid on;



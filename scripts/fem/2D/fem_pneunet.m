clr;
%% settings
Mat = Dragonskin10();
P0  =  7.5*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;
MinH  = 2;
MaxH  = 3;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh = Blender(msh,'Rotate',{'x',90});
msh = msh.show();
pause(2);

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/35,'Linestyle','none','StressNorm',1e-3);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,5]),[1,1]);

id = fem.FindEdges('AllHole');
fem = fem.AddConstraint('Pressure',id,[P0,0]);

%% add output nodes
id = fem.FindNodes('Right');
fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = Mat;

%% solve
fem.solve();

%% post-processing
t = fem.Log{1,2};
Ux = fem.Log{2,2};
Uy = fem.Log{3,2};
Ve = fem.Log{12,2};
N0 = fem.get('Node0');

figure(103); cla; subplot(1,2,1);
for ii = 1:1:size(Ux,2)
    Nx = N0(id,1) + Ux(:,ii);
    Ny = N0(id,2) + Uy(:,ii);
    plot(Nx,Ny,'Linewidth',2,...
        'Color',col(4,ii/(1.05*size(Ux,2))));
    hold on;
    Y{ii} = 1e-3*[Nx,Ny];
end

axis([-20 100 -20 130]); axis equal;
xaxis('$x$-dimension','mm');
yaxis('$y$-dimension','mm');
set(gca,'linewidth',1)
grid on;

subplot(1,2,2);
plot([0;t],[0;Ve],'Linewidth',2,...
        'Color',col(1));
xaxis('Normalized loading','-');
yaxis('Potential energy','J');
set(gca,'linewidth',1)
grid on;

%% shape fitting

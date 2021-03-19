clr;
%%
I = imread('bitmap.png');
Ir = imcomplement(im2gray(I));
Ir(:,1:5) = [];
msh = GenerateMeshImage(Ir,[-28 28 0 62],0.08,[1,1,40]);

msh.show();

%%
fem = Fem(msh,'Nonlinear',true,'TimeStep',1/15,'FilterRadius',60/15,...
    'Movie',0,'Linestyle','none','Linestyle0','none');

%% assign boundary conditions to reduced fem
id = fem.FindEdges('BoxHole',[34,45,3,60]);
fem = fem.AddConstraint('Pressure',id,[75*kpa,0]);

id = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Box',[27,29,0,62]);
fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material to reduced fem
fem.Material = TPU90;

%% solve final finite-element problem
fem.solve();

%% post-processing
Ux = fem.Log{1,2};
Uy = fem.Log{2,2};
N0 = fem.get('Node0');

figure(103); cla;
for ii = 15:1:size(Ux,2)
    Nx = N0(id,1) + Ux(:,ii);
    Ny = N0(id,2) + Uy(:,ii);
    plot(Nx,Ny,'.','Linewidth',1.5,...
        'Color',col(1,ii/size(Ux,2)));
    hold on;
end
axis equal;

%% fitting
ModelFit([Nx-0.028e3,Ny-1.5]*1e-3,0.062,[0,1,0,1,0,0],4);
%BezierFit([Nx,Ny],{1:6,101:110});
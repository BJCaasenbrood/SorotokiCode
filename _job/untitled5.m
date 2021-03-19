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
id = fem.FindEdges('Hole',[4,7,10,13,16,19,22,25]);

fem = fem.AddConstraint('Pressure',id,[55*kpa,0]);

id = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Support',id,[1,1]);

%% assign material to reduced fem
fem.Material = TPU90;

%% solve final finite-element problem
fem.solve();
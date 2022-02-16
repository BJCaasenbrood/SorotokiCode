clr;

%% generate mesh
Simp = 0.09;
W = 10*4.35;
H = 10*14.35; 

msh = Mesh('SWR_V3.png','BdBox',[0,W,0,H],'SimplifyTol',Simp,...
    'Hmesh',[1.1,1.5,7]);

msh = msh.generate();
msh.show(); 

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/10,'Linestyle','none',...
    'MovieAxis',[-1 45 -1  170]);

%% assign material
fem.Material = NeoHookeanMaterial(60,0.3);

%% add boundary constraint
id = fem.FindEdges('Hole',2);
fem = fem.AddConstraint('Support',id{:},[1,1]);

id = fem.FindEdges('Hole',12);
fem = fem.AddConstraint('Support',id{:},[1,0]);

id = fem.FindEdges('Hole',12);
fem = fem.AddConstraint('Spring',id{:},[0,0.01]);

fem = fem.AddConstraint('Output',id{:},[1,0]);

P0 = -25*kpa;
% 
id = fem.FindEdges('Hole',2:11);
fem = fem.AddConstraint('Pressure',id,P0);

%% solve
f = figure(101);
fem.solve();

%%

X = fem.Log.t*P0/kpa;
Y = mean(0.2*fem.Log.Uy,1) + 0.00*(-0.5*X).^2;

X(end + 1) = X(end) -2.5;
X(end + 1) = X(end) -2.5;
X(end + 1) = X(end) -2.5;
X(end + 1) = X(end) -2.5;
dY = mean(diff(Y));

Y(end + 1) = Y(end) + dY + 0*1.75;
Y(end + 1) = Y(end) + dY + 0*.3;
Y(end + 1) = Y(end) + dY + 0;
Y(end + 1) = Y(end) + dY + 0;

%figure(103);
plot(-X,Y,'-s','LineW',3);


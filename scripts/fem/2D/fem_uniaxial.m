clr;
%% simulation settings
H  = 20;       % height of specimen
W  = 20;       % width of specimen
dL = H*5;      % elongation of specimen
 
%% signed distance function (SDF)
sdf = sRectangle(0,W,0,H);
 
%% generate mesh
msh = Mesh(sdf,'Quads',[2,2]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/20,'PrescribedDisplacement',true,...
              'Linestyle','none','ColorAxis',[0 1],'SolverPlot',true);

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,dL]);

fem = fem.AddConstraint('Output',fem.FindNodes('Location',[0,H]),[0,0]);
%% assign material
%fem.Material = Ecoflex0050;
fem.Material = Ecoflex0030_Ogden;

%% solving extension
fem.solve();

%% plotting
figure(102)
lam = (20+fem.Log.Uy)/H;
svm = fem.Log.Svm;

hold on;
plot(lam,svm,'-o','linewidth',2);
xaxis('Elongation strain $\lambda$','-')
yaxis('Nominal stress','Mpa')
clr;
%% parameters
Elongation = 500;

%% generate mesh from sdf
F = [1,2,3,4];
V = [0,0;5,0;5,5;0,5];

msh = Mesh(V,F);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'PrescribedDisplacement',true,...
              'Linestyle','none','ColorAxis',[0 1]);

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[5*Elongation/100,0]);
fem = fem.AddConstraint('Output',fem.FindNodes('Location',[5,0]),[0,0]);
%% assign material
fem.Material = Ecoflex0050(0.1);

%% solving extension
fem.solve();
lam = fem.Log{1,2}*Elongation;
svm = fem.Log{8,2};

%% plot
figure(102)

hold on;
plot(lam,svm,'linewidth',2);
xaxis('Elongation strain $\lambda$','-')
yaxis('Nominal stress','Mpa')
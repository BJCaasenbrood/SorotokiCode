%clr;
%% parameters
Stretch = [-0.75,3];

%% generate mesh
F = [1,2,3,4];
V = [0,0;5,0;5,5;0,5];

msh = Mesh(V,F);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'PrescribedDisplacement',true,...
              'Linestyle','none','ColorAxis',[0 1]);

%% add boundary conditions
fem = fem.addSupport(fem.FindNodes('Left'),[1,0]);
fem = fem.addSupport(fem.FindNodes('SW'),[1,1]);
fem = fem.addLoad(fem.FindNodes('Right'),[msh.BdBox(2)*Stretch(2),0]);
%fem = fem.addOutput(fem.FindNodes('Location',[5,0]),[0,0]);

%% assign material
fem.Material = Ecoflex0050();

%% solving extension
fem.solve();
t1 = fem.Log{1,2}*Stretch(2);
s1 = fem.Log{8,2};

%% solving compression
fem = fem.set('TimeStep',1/20);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[msh.BdBox(2)*Stretch(1),0]);

fem.solve();
t2 = fem.Log{1,2}*Stretch(1);
s2 = -fem.Log{8,2};

%% stitch data
lam = [flipud(t2);0;t1];
svm = [fliplr(s2),0,s1];

%% plot
figure(102)

hold on;
plot(lam,svm,'linewidth',2);
xaxis('Elongation strain $\lambda$','-')
yaxis('Nominal stress','Mpa')
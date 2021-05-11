%clr;
%% parameters
Stretch = [-0.75,2];

%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,5);

msh = Mesh(sdf,'BdBox',[0,5,0,5],'Quads',[2,2]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'PrescribedDisplacement',true,...
              'Linestyle','none','ColorAxis',[0 1]);

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[msh.BdBox(2)*Stretch(2),0]);
fem = fem.AddConstraint('Output',fem.FindNodes('Location',[5,0]),[0,0]);
%% assign material
fem.Material = Ecoflex0050;

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
lam = [fliplr(t2),0,t1];
svm = [fliplr(s2),0,s1];

%% plot
figure(102)

hold on;
plot(lam,svm,'linewidth',2);
xaxis('Elongation strain $\lambda$','-')
yaxis('Nominal stress','Mpa')
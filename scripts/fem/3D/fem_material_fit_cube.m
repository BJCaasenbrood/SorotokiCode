clr;
%% generate mesh from STL
msh = Mesh('Cube0.stl','Hmesh',[1,1,1]);
msh = msh.show();

%% parameters
Stretch = [-1.5,6];
Material = YeohMaterial('C1',0.3,'C2',0.01);

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/30,'PrescribedDisplacement',true,...
    'Linestyle','-','ColorAxis',[-2,2]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Front'),[0,1,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,0,msh.BdBox(2)*Stretch(2)]);
fem = fem.AddConstraint('Output',fem.FindNodes('Location',[-0.5,-0.5,0]),[0,0,0]);

%% select material
fem.Material = Material;

%% solving
fem.solve();
t1 = fem.Log{1,2}*Stretch(2);
s1 = fem.Log{8,2};

%% solving compression
fem = fem.set('TimeStep',1/20);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,0,msh.BdBox(2)*Stretch(1)]);

fem.solve();
t2 = fem.Log{1,2}*Stretch(1);
s2 = -fem.Log{8,2};

%% stitch data
lam = [flipud(t2(:));0;t1(:)];
svm = [flipud(s2(:));0;s1(:)];

%% plot
figure(102); overwrite_colors;
hold on;
plot(lam,svm,'linewidth',2);
xaxis('Elongation strain $\lambda$','-');
yaxis('Nominal stress','Mpa');
grid on; box on; 
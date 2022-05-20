%% simulation settings
H  = 20;       % height of specimen
W  = 20;       % width of specimen
dL = H*5;      % elongation of specimen

%% signed distance function (SDF)
sdf = sRectangle(0,W,0,H);

%% generate mesh
msh = Mesh(sdf,'Quads',[20,20]);
msh = msh.generate();

%% show SDF and mesh
figure(101);
subplot(1,2,1); sdf.show();
subplot(1,2,2); msh.show();

%% generate fem model
clf;
fem = Fem(msh,'TimeStep',1/25,'Linestyle','none','ColorAxis',[0,1]);

%% adding boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Top'),[0,dL]);

%% outputs nodal data in fem.Log
fem = fem.AddConstraint('Output',fem.FindNodes('NW'),[0,0]);

%% adding material
fem.Material = Ecoflex0030;

%% solving
fem = fem.solve();

%% post-processing data and plotting
Exx = (dL/H)*fem.Log.t;
Svm = fem.Log.Svm;

figure(102);
subplot(1,2,1); fem.show('Uy');
subplot(1,2,2); plot(Exx,Svm,'-o','linewidth',2)
xlabel('Uni-axial strain (-)','interpreter','latex','fontsize',20);
ylabel('Von Mises stress (MPa)','interpreter','latex','fontsize',20);
set(gca,'linewidth',1.5); grid on; axis tight;

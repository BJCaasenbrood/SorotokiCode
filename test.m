clr;
H  = 20;       % height of specimen
W  = 20;       % (half) width of specimen
dL = H*4;      % elongation of specimen
sdf = sRectangle(0,W,0,H);
msh = Mesh(sdf,'Quads',[5,5]);
msh = msh.generate();

figure(101);
subplot(1,2,1); sdf.show();
subplot(1,2,2); msh.show();
fem = Fem(msh,'TimeStep',1/15,'PrescribedDisplacement',true,...
'Linestyle','none','SolverPlotType','Uy');
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,dL]);
fem = fem.AddConstraint('Output',fem.FindNodes('NW'),[0,0]);
fem.Material = Ecoflex0030(0.1);
fem = fem.solve();
%% post-processing data and plotting
Eyy = (dL/H)*fem.Log.t + 1;
Syy = fem.Log.Svm;

% exact S_xx solution
Eyy_exact = linspace(1,1+(dL/H),500);
Syy_exact = @(x) 2*(x.^2 - 1./x).*fem.Material.dWdI(x.^2 + 2./x);

figure(102);
subplot(1,2,1); fem.show('Uy'); 
subplot(1,2,2); plot(Eyy,Syy,'o-','Color',col(1),'linewidth',2); 
subplot(1,2,2); hold on; plot(Eyy_exact,Syy_exact(Eyy_exact),...
    '--','Color',col(1),'linewidth',2);

xlabel('Uni-axial strain (-)','interpreter','latex','fontsize',20);
ylabel('Von Mises stress (MPa)','interpreter','latex','fontsize',20);
legend('FEM','Exact','Location','Northwest','fontsize',18);
grid on; axis tight;


clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,20,0,2);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,2],'NElem',150);
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/200,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('NE'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[-4,1e-3]);

%% add logger nodes
fem = fem.AddConstraint('Output',fem.FindNodes('SE'),[0,0]);

%% assign material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

%% plot force-displacement relation
figure(101);
subplot(2,1,1); fem.show();
subplot(2,1,2); plot(fem.Log{3},fem.Log{6},'linewidth',2,'Color',col(2));
xlabel('Displacement (mm)','interpreter','latex','fontsize',12);
ylabel('Reaction force (N)','interpreter','latex','fontsize',12);
grid on; set(gca,'linewidth',1);

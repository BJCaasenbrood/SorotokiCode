clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,20,0,1);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,1],'NElem',150);
msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/200,'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[-2,0]);

%% add logger nodes
fem = fem.AddConstraint('Output',fem.FindNodes('SE'),[0,0]);

%% assign material
fem.Material = Dragonskin10A;

%% solving
fem.solve();

%% plot force-displacement relation
figure(101);
subplot(2,1,1); fem.show();
subplot(2,1,2); plot(fem.Log{3,2},fem.Log{6,2}*1e3,'linewidth',2,'Color',col(2));
xlabel('Displacement (mm)','interpreter','latex','fontsize',12);
ylabel('Reaction force (mN)','interpreter','latex','fontsize',12);
grid on; set(gca,'linewidth',1);

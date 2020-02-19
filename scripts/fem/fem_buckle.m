clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,20,0,2);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,2],...
              'NElem',50,...
              'MaxIteration',150,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/200,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('NE'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[-10,0]);

%% add logger nodes
fem = fem.AddConstraint('Output',fem.FindNodes('SE'),[0,0]);

fem.Material = Dragonskin10A;

%% solving
fem.solve();

%% plot force-displacement relation
figure(102);
plot(fem.Log{3},fem.Log{6}*1e3,'linewidth',2,'Color',col(1));
axis tight; grid on;
xlabel('Displacement (mm)','interpreter','latex');
ylabel('Reaction force (mN)','interpreter','latex');
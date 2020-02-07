%% set signed distance function
sdf = @(x) dRectangle(x,0,8,0,2);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,8,0,2],...
              'NElem',500,...
              'MaxIteration',150,...
              'ShowMeshing',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'FilterRadius',0.5,...
              'VolumeInfill',0.3,...
              'Penal',4,...
              'OptimizationProblem','Compliance',...
              'PrescribedDisplacement',false,...
              'Nonlinear',true);

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[4,2],3); 
fem = fem.AddConstraint('Load',id,[0,-2e-4]);

%% material
fem.Material = Ecoflex0030;
% fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
%     'D1',1.5,'D2',2.0,'D3',1.0);

%% set density
fem = fem.initialTopology('Equidistance',[3,1],.4);

%% solving
fem.optimize();

%% former
fem.former();
fem.showTopo(0.15);

% fem = fem.AddConstraint('Load',id,[0,-3e-3]);
% fem = fem.set('TimeStep',1/15);
% fem.solve();
% fem.show('Sxy')

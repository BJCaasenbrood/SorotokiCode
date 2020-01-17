%% set signed distance function
sdf = @(x) dRectangle(x,0,6,0,2);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,6,0,2],...
              'NElem',1000,...
              'MaxIteration',150,...
              'ShowMeshing',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'FilterRadius',0.3,...
              'VolumeInfill',0.3,...
              'Penal',1,...
              'OptimizationProblem','Compliance',...
              'PrescribedDisplacement',false,...
              'ReflectionPlane',[-1,0],...
              'Nonlinear',false);

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,0]);
id = fem.FindNodes('SE'); 
fem = fem.AddConstraint('Support',id,[0,1]);

id = fem.FindNodes('NW',[4,2],1); 
fem = fem.AddConstraint('Load',id,[0,-1e-4]);

%% material
fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',1.5,'D2',2.0,'D3',1.0);

%% set density
fem = fem.initialTopology('Equidistance',[2,1],.4);

%% solving
fem.optimize();

%% former
fem.former();
fem.showTopo(0.15);
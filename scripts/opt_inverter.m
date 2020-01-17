%% set signed distance function
sdf = @(x) dRectangle(x,0,5,0,2);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,2],...
              'NElem',500,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/10,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.3,...
              'Penal',1,...
              'PenalMax',4,...
              'PrescribedDisplacement',true,...
              'FilterRadius',0.2,...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Top'); 
fem = fem.AddConstraint('Support',id,[0,1]);

id = fem.FindNodes('Location',[0,0],5);
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[0,2],1);
fem = fem.AddConstraint('Load',id,[1,0]);
fem = fem.AddConstraint('Spring',id,[.1,0]);

id = fem.FindNodes('Location',[5,2],1);
fem = fem.AddConstraint('Output',id,[-1,0]);
fem = fem.AddConstraint('Spring',id,[.1,0]);

%% set density
fem = fem.initialTopology([1,1],.5);

%% material
% fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
%     'D1',15,'D2',20,'D3',10);

% fem.Material = YeohMaterial('C1',1,'C2',0,'C3',0,...
%     'D1',1.0,'D2',1.0,'D3',1.0);

fem.Material = LinearMaterial('E',10,'Nu',0.49);

% fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
%     'D1',1.5,'D2',2.0,'D3',1.0);


fem.show('E');

%% solving
fem.optimize();

%% former
fem.former();
fem.showTopo(0.15);



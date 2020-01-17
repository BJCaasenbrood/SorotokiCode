%% set signed distance function
Height = 12;
Width = 4;
Depth = 5;

sdf = @(x) Hanger(x,Height,Width,Depth);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,Width+Depth,-Height,0],...
              'NElem',750,...
              'MaxIteration',150,...
              'ShowMeshing',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'FilterRadius',0.75,...
              'VolumeInfill',0.5,...
              'Penal',1,...
              'OptimizationProblem','Compliance',...
              'PrescribedDisplacement',false,...
              'Nonlinear',false);

%% add constraint
id = fem.FindNodes('Top'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Location',[Depth+Width,Width/2 - Height],1); 
fem = fem.AddConstraint('Load',id,[0,-1e-4]);

%% material
fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
     'D1',1.5,'D2',2.0,'D3',1.0);

%fem.Material = LinearMaterial('E',1,'Nu',0.49);

%% solving
fem.optimize();

function Dist = Hanger(X,H,W,D)
R1 = dRectangle(X,0,W,-H,0);
R2 = dRectangle(X,W,W+D,-H,W-H);
Dist = dUnion(R1,R2);
end

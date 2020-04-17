clr;
%% set signed distance function
sdf = @(x) SDF(x);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,30,0,50],'NElem',1e3);    
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'TimeStep',1/3,...
              'ResidualNorm',1e-3,...
              'FilterRadius',0.3,...
              'VolumeInfill',0.3,...
              'ChangeMax',Inf,...
              'Penal',4,...
              'Nonlinear',false,...
              'OptimizationProblem','Compliance');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Top'); 
fem = fem.AddConstraint('Load',id,[0,-1e-3]);

%% material
fem.Material = Dragonskin10A;

%% set density
fem.Density = ones(fem.NElem,1);

%% solving
fem.optimize();

function D = SDF(x)

R1 = dRectangle(x,0,0,2,1);
C1 = dCircle(x,0,0.5,0.1);

D = R1;

end

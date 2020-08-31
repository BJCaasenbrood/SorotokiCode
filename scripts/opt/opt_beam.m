clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,4,0,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,4,0,2],'Quads',500);      
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'TimeStep',1/3,'ResidualNorm',1e-3,'FilterRadius',0.3,...
        'VolumeInfill',0.3,'ChangeMax',Inf,'Penal',4,'Nonlinear',true,...
              'OptimizationProblem','Compliance');

%% set symmetry          
fem = fem.set('Periodic',[1, 0],'ReflectionPlane',[-1,0]);          
          
%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Support',id,[1,0]);

id = fem.FindNodes('Location',[4,2]); 
fem = fem.AddConstraint('Load',id,[0,-2e-5]);

%% material
fem.Material = Ecoflex0030;

%% set density
fem.Density = ones(fem.NElem,1);

%% solving
fem.optimize();
fem.show('ISO',0.3);

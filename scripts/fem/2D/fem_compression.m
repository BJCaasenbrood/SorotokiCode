clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,10],'NElem',500,'MaxIteration',800);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'Nonlinear',true,...
              'PrescribedDisplacement',true,'Penal',2,...
               'Movie',true,'MovieAxis',[-0.1 10 0 13]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Line',[0 4 10 10]),[1,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Line',[0 10 10 10]),[0,2]);

%% assign material
fem.Material = Ecoflex0030;

Pc = fem.Mesh.get('Center');
fem.Density = ones(500,1);
fem.Density(Pc(:,2) < 5) = 0.5;

%% solving
figure(101);
fem.solve();
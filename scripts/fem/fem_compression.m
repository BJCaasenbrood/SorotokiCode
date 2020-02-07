%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,10],...
              'NElem',250,...
              'MaxIteration',500,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/25,...
              'ResidualNorm',1e-3,...
              'DisplaceNorm',1e-3,...
              'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Line',[0 4 10 10]),[0,-4]);
fem = fem.AddConstraint('Support',fem.FindNodes('Line',[0 4 10 10]),[1,0]);

%fem.Material = Ecoflex0030('Yeoh');
fem.Material = NeoHookeanMaterial('E',5,'Nu',0.25);

%% solving
fem.solve();

%% plotting
figure(101); clf;
fem.show('Svm');
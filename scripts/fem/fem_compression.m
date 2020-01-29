%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,10],...
              'NElem',150,...
              'MaxIteration',500,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
msh.show();

fem = Fem(msh);
fem = fem.set('TimeStep',1/15,...
              'ResidualNorm',1e-6,...
              'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Line',[0 2 10 10]),[0,-3]);

fem.Material = Ecoflex0030('Yeoh');
%fem.Material = NeoHookeanMaterial('E',5,'Nu',0.45);

%% solving
fem.solve();

%% plotting
figure(101); clf;
fem.show('Svm');
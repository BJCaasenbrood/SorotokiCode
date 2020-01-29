clc; clear; close all;
%% set signed distance function
sdf = @(x) dRectangle(x,0,1,0,.5);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,1,0,.5],...
              'NElem',500,...
              'MaxIteration',150,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/100,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'PrescribedDisplacement',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[5,0]);

% fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
%      'D1',1.5,'D2',1,'D3',1.0);
%fem.Material = MooneyMaterial('C10',2,'K',300);

fem.Material = Ecoflex0030;

%% solving
fem.solve();

%% plotting
figure(101); clf;
fem.show('Svm');
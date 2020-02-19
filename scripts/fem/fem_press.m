close all; clc; clear;
%% set signed distance function
sdf = @(x) dCircle(x,0,0,1);
%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[-1,1,-1,1],...
              'Center',SeedList);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/15,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'FilterRadius',1e-7,...
              'Penal',4,...
              'VolumetricPressure',true,...
              'PrescribedDisplacement',false);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom',2),[1,1]);
fem = fem.AddConstraint('PressureCell',fem.FindElements(...
    'Location',[0,0],1),[1e-4,0]);

fem.Material = Dragonskin10A;

%% generate void region
%fem = fem.initialTopology('Hole',[0,0],.75);
% 
id = fem.FindElements('Location',[0,0]);
fem.Density = ones(fem.NElem,1);
fem.Density(id) = 1e-2;

%% solving
fem.solve();


function P = SeedList

P =[0.5735, -0.4550;
   -0.7151,  0.1562;
   -0.3240,  0.6566;
    0.0027, -0.7319;
    0.0005,  0.0007;
    0.3129,  0.6621;
    0.7132,  0.1658;
   -0.5691, -0.4601];

end

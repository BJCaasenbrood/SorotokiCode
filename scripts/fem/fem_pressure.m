close all; clc; clear;
%% set signed distance function
sdf = @(x) Square(x,20,10);
%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,10],...
              'NElem',500,...
              'MaxIteration',300,...
              'ShowMeshing',false,...
              'Triangulate',false);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/15,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'Penal',4,...
              'VolumetricPressure',true,...
              'PrescribedDisplacement',false);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[1,1]);
fem = fem.AddConstraint('PressureCell',fem.FindElements(...
    'Location',[5,5],1),[4e-4,0]);

fem.Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
      'D1',10,'D2',10,'D3',10);

%% generate void region
f = @(x) dRectangle(x,2,18,2,8);
d = f(fem.Center); id = find(d(:,end) < 0);
fem.Density(id) = 1e-3;

%% solving
fem.solve();

function D = Square(P,W,H)
D = dRectangle(P,0,W,0,H);
end
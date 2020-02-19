close all; clc; clear;
%% set signed distance function
sdf = @(x) SDF(x,20);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,20],...
              'NElem',500);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/10,...
              'ResidualNorm',1e-3,...
              'FilterRadius',1,...
              'Penal',4,...
              'VolumetricPressure',true,...
              'Nonlinear',false);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('PressureCell',fem.FindElements(...
    'Location',[5,5],1),[1e-3,0]);

fem.Material = Dragonskin10A;

%% generate void region
f = @(x) dCircle(x,0,0,18);
d = f(fem.Center); id = find(d(:,end) < 0);
fem.Density(id) = 1e-2;

%% solving
fem.solve();

function Dist = SDF(P,R)
R1 = dRectangle(P,0,R,0,R);
C1 = dCircle(P,0,0,R);
Dist = dIntersect(C1,R1);
end

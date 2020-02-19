clr;
%% set signed distance function
sdf = @(x) SDF(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,20],...
              'NElem',300);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/2,...
              'ResidualNorm',1e-3,...
              'VolumeInfill',0.3,...
              'VolumetricPressure',true,...
              'FilterRadius',4,...
              'Type','PlaneStrain',...
              'Penal',4,...
              'Nonlinear',false,...
              'MaxIterationMMA',75,...
              'Movie',false,...
              'MovieAxis',[0 100 0 20],...
              'ReflectionPlane',[-1 0],...
              'Repeat',[1,1,1],...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);
%id = fem.FindNodes('Bottom'); 
%fem = fem.AddConstraint('Support',id,[0,1]);
%fem = fem.AddConstraint('Spring',id,[0,1e-3]);
% id = fem.FindNodes('Line',[0 1 0 0]); 
% fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Output',id,[1,-0.1]);
fem = fem.AddConstraint('Spring',id,[1e-3,1e-3]);

% id = fem.FindElements('Bottom'); 
% fem = fem.AddConstraint('Spring',id,[1e-5,1e-5]);

id = fem.FindElements('Location',[5,10],1);
fem = fem.AddConstraint('PressureCell',id,[0.1,0]);

%% set density
fem = fem.initialTopology('Hole',[5,10],1);

%% material
fem.Material = Dragonskin10A;

%fem.Material = LinearMaterial('E',3,'Nu',0.49);


%% solving
fem.optimize();

%% show interpolated surface
fem.set('Repeat',[1,1,1],'ReflectionPlane',[1,1]);
fem.show('ISO');

function Dist = SDF(P)
  R1 = dRectangle(P,0,10,0,20);
  R3 = dRectangle(P,9,12,4,22);
  C1 = dCircle(P,10,4.5,1);
  Dist = dDiff(dDiff(R1,R3),C1);
end


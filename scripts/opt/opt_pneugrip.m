clr;
%% set signed distance function
sdf = @(x) PneuGrip(x);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-5,13,0,5],'NElem',3e3);%'Quads',[50 25]);
msh = msh.generate();

%% generate fem model
fem = Fem(msh,'VolumeInfill',0.25,'Penal',1,'FilterRadius',0.5,...
              'Nonlinear',false,'TimeStep',1/3,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',125,'ChangeMax',0.02,'Movie',false);
          
%% set spatial settings
fem = fem.set('ReflectionPlane',[0 -1]);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,1]);

id = fem.FindNodes('Location',[13,4.5],1);
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindElements('Location',[-5,5],1);
fem = fem.AddConstraint('PressureCell',id,[-1*kpa,0]);

%% set density
fem = fem.initialTopology('Hole',[-5,5],2);

%% material
fem.Material = Ecoflex0050(0.75);

%% solving
fem.optimize();
fem.show('ISO',0.15);

%% convert topology result to mesh
mshr = fem.exportMesh(0.2,0.035,[1.0,0.3,5]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/25,'FilterRadius',2/15,...
    'MovieAxis',[-75 170 -140 40],'Movie',0,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Left'); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindEdges('EdgeSelect',[10.7,4.44],60);
femr = femr.AddConstraint('Pressure',id,[4*kpa,0]);

id = femr.FindNodes('Location',[17.64,3.05],1);
femr = femr.AddConstraint('Support',id,[1,0]);

id = femr.FindNodes('Location',[17.64,6.762],1);
femr = femr.AddConstraint('Support',id,[1,0]);

id = femr.FindNodes('Location',[17.64,6.762],1);
femr = femr.AddConstraint('Support',id,[1,0]);

id = femr.FindNodes('Location',[0.2352,8.705],1);
femr = femr.AddConstraint('Support',id,[1,1]);
id = femr.FindNodes('Location',[0.2352,1.11],1);
femr = femr.AddConstraint('Support',id,[1,1]);
%% assign material to reduced fem
D = 25; % compress. factor (more stable)
femr.Material = Ecoflex0050(D);

%% solve final finite-element problem
femr.solve();

function Dist = PneuGrip(P)
  R1 = dRectangle(P,-5,13,0,5);
  R2 = dRectangle(P,10,15,3,7);
  R3 = dRectangle(P,10,15,0,1);
  C2 = dCircle(P,10,5,2);
  C3 = dCircle(P,10,0,1);
  Dist = dDiff(dDiff(dDiff(dDiff(R1,C2),R2),R3),C3);
end


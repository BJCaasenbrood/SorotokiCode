clr;
%% generate mesh from sdf
W = 20;  % width cell
H = 40;  % width cell
D = 2;   % inter distance

sdf = @(x) PneuNet(x,W,H,D,W);

msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',500);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.4,'Penal',4,'FilterRadius',H/20,...
              'Nonlinear',1,'TimeStep',1/30,'ReflectionPlane',[0,0],...
              'OptimizationProblem','Compliant','Linestyle','None',...
              'MaxIterationMMA',70,'ChangeMax',0.1,'Movie',0);

%% set spatial settings
fem = fem.set('Periodic',[1/2, 0],'Repeat',[ones(9,1)]);

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id  = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);
id  = fem.FindElements('Location',[W/2,0.5*H],1);
fem = fem.AddConstraint('PressureCell',id,10*kpa);

%% set density
fem = fem.initialTopology('Hole',[W/2,0.5*H],0.85);

%% material
fem.Material = Ecoflex0030(2);

%% solving
fem.optimize();
fem.show('ISO',0.2);

%% convert topology result to mesh
ISO  = 0.3;
Simp = 0.05;
GrowH = 1.5;
MinH = 2.5;
MaxH = 40;

mshr = fem.exportMesh(ISO,Simp,[GrowH,MinH,MaxH]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/120,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Left'); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindEdges('TopHole');
femr = femr.AddConstraint('Pressure',id,15*kpa);

id = femr.FindNodes('Bottom');
femr = femr.AddConstraint('Output',id,[0,0]);

%% assign material to reduced fem
D = 25; % compress. factor (more stable)
femr.Material = Ecoflex0030(D);

%% solve final finite-element problem
femr.solve();

%% post-process curvature data
t = femr.Log.t; 
close all;
figure(105);

for ii = 1:fps(t,60):numel(t)
    N = femr.Log.Node{ii};

    femr.set('Node',N);
    femr.show('Field',femr.Log.Stress{ii});
    axis([-70 270 -200 50]);
    background(gitpage);
    
    drawnow();
end

function Dist = PneuNet(P,W,H,E,T)
R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,-W/2,E,T,H+H/2);
R3 = dRectangle(P,W-E,W+W/2,T,H+H/2);
C1 = dCircle(P,0,T + 0.5,1);
C2 = dCircle(P,W,T + 0.5,1);
Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end
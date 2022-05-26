clr;
%% generate mesh from sdf
W = 30;  % width cell
H = 55;  % width cell
D = 4;   % inter distance

sdf = @(x) PneuNet(x,W,H,D,W);

msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',700);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',2,'FilterRadius',H/15,...
              'Nonlinear',false,'TimeStep',1/3,'ReflectionPlane',[0,1],...
              'OptimizationProblem','Compliant','Linestyle','None',...
              'MaxIterationMMA',50,'ChangeMax',0.03,'Movie',0);

%% set spatial settings
fem = fem.set('Periodic',[1/2, 0],'Repeat',[ones(1,1)]);

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id  = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);
id  = fem.FindElements('Location',[W/2,0.5*H],1);
fem = fem.AddConstraint('PressureCell',id,3*kpa);

%% set density
fem = fem.initialTopology('Hole',[W/2,0.5*H],0.85);

%% material
fem.Material = Ecoflex0030();

%% solving
fem.optimize();
fem.show('ISO',0.25);

%% convert topology result to mesh
ISO  = 0.3;
Simp = 0.05;
GrowH = 1;
MinH = 4;
MaxH = 30;

mshr = fem.exportMesh(ISO,Simp,[GrowH,MinH,MaxH]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/50,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Left'); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindEdges('BoxHole',[0,200,60,100]);
femr = femr.AddConstraint('Pressure',id,6*kpa);

id = femr.FindEdges('BoxHole',[200,450,0,50]);
femr = femr.AddConstraint('Pressure',id,6*kpa);

id = femr.FindNodes('Bottom');
femr = femr.AddConstraint('Output',id,[0,0]);

%% assign material to reduced fem
D = 20; % compress. factor (more stable)
femr.Material = Ecoflex0030(D);

%% solve final finite-element problem
femr.solve();

%% movie
t = femr.Log.t; 
close all;
figure(105);

for ii = 1:fps(t,60):numel(t)
    N = femr.Log.Node{ii};

    femr.set('Node',N);
    femr.show('Field',femr.Log.Stress{ii});
    axis([-20 350 -200 120]);
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
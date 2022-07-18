clr;
%% generate mesh from sdf
W = 20;  % width cell
H = 55;  % width cell
D = 3;   % inter distance

sdf = @(x) PneuNet(x,W,H,D,W);

msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',3e3);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.33,'Penal',2,'FilterRadius',H/20,...
              'Nonlinear',0,'TimeStep',1/10,'ReflectionPlane',[0,1],...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',50,'ChangeMax',0.05);

%% set spatial settings
fem = fem.set('Periodic',[1/2, 0],'Repeat',[ones(10,1)]);

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id  = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);
id  = fem.FindElements('Location',[W/2,0.3*H],1);
fem = fem.AddConstraint('PressureCell',id,0.1*kpa);

%% set density
fem = fem.initialTopology('Hole',[W/2,0.3*H],.5);

%% material
fem.Material = Dragonskin10();

%% solving
fem.optimize();
fem.show('ISO',0.3); drawnow;

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

id = femr.FindEdges('TopHole');
femr = femr.AddConstraint('Pressure',id,20*kpa);

%% assign material to reduced fem
femr.Material = fem.Material;

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
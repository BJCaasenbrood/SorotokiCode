clr;
%% generate mesh from sdf
W = 30;  % width cell
H = 75;  % width cell
D = 2;   % inter distance

sdf = @(x) PneuNet(x,W,H,D,W);

msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',800);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.35,'Penal',4,'FilterRadius',H/15,...
              'Nonlinear',false,'TimeStep',1/3,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',25,'ChangeMax',0.05,'Movie',0);

%% set spatial settings
fem = fem.set('Periodic',[1/2, 0],'Repeat',ones(7,1));

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);
id = fem.FindElements('Location',[W/2,0.625*H],1);
fem = fem.AddConstraint('PressureCell',id,[5*kpa,0]);

%% set density
fem = fem.initialTopology('Hole',[W/2,0.625*H],0.85);

%% material
fem.Material = Ecoflex0030();

%% solving
fem.optimize();
fem.show('ISO',0.25);

%% convert topology result to mesh
ISO  = 0.3;
Simp = 0.05;
GrowH = 1.5;
MinH = 2.5;
MaxH = 40;

mshr = fem.exportMesh(ISO,Simp,[GrowH,MinH,MaxH]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/35,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Left'); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindEdges('TopHole');
femr = femr.AddConstraint('Pressure',id,[10*kpa,0]);

id = femr.FindNodes('Bottom');
femr = femr.AddConstraint('Output',id,[0,0]);

%% assign material to reduced fem
D = 25; % compress. factor (more stable)
femr.Material = Ecoflex0030(15);

%% solve final finite-element problem
femr.solve();

%% post-process curvature data
Ux = femr.Log{2,2};
Uy = femr.Log{3,2};
N0 = femr.get('Node0');

figure(103); cla;
for ii = 1:1:size(Ux,2)
    Nx = N0(id,1) + Ux(:,ii);
    Ny = N0(id,2) + Uy(:,ii);
    plot(Nx,Ny,'Linewidth',2,...
        'Color',col(4,ii/(1.05*size(Ux,2))));
    hold on;
end

function Dist = PneuNet(P,W,H,E,T)
R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,-W/2,E,T,H+H/2);
R3 = dRectangle(P,W-E,W+W/2,T,H+H/2);
C1 = dCircle(P,0,T + 0.5,1);
C2 = dCircle(P,W,T + 0.5,1);
Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end
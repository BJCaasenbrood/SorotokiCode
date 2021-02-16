clr;
%% generate mesh from sdf
W = 20;  % width cell
H = 45;  % width cell
D = 2;   % inter distance

sdf = @(x) PneuNet(x,W,H,D,W);

msh = Mesh(sdf,'BdBox',[0,W,0,H],'Quads',[25 50]);
msh = msh.generate();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',H/15,...
              'Nonlinear',false,'TimeStep',1/3,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',15,'ChangeMax',0.15,'Movie',0);

%% set spatial settings
fem = fem.set('Periodic',[1/2, 0],'Repeat',ones(7,1));

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);

id = fem.FindElements('Location',[W/2,0.625*H],1);
fem = fem.AddConstraint('PressureCell',id,[4*kpa,0]);

%% set density
fem = fem.initialTopology('Hole',[W/2,0.625*H],0.15);

%% material
fem.Material = Ecoflex0030(0.75);

%% solving
fem.optimize();
fem.show('ISO',0.3);

%% convert topology result to mesh
mshr = fem.exportMesh(0.3,0.05,[1.1,1.85,40]); 
mshr.show(); pause(2);

femr = Fem(mshr,'Nonlinear',true,'TimeStep',1/25,'FilterRadius',H/15,...
    'Movie',0,'Linestyle','none');

%% assign boundary conditions to reduced fem
id = femr.FindNodes('Left'); 
femr = femr.AddConstraint('Support',id,[1,1]);

id = femr.FindEdges('TopHole');
femr = femr.AddConstraint('Pressure',id,[4*kpa,0]);

id = femr.FindNodes('Bottom');
femr = femr.AddConstraint('Output',id,[0,0]);

%% assign material to reduced fem
D = 25; % compress. factor (more stable)
femr.Material = Ecoflex0030(D);

%% solve final finite-element problem
femr.solve();

%% post-process curvature data
Ux = femr.Log{1,2};
Uy = femr.Log{2,2};
N0 = femr.get('Node0');

figure(103); cla;
for ii = 1:1:size(Ux,2)
    Nx = N0(id,1) + Ux(:,ii);
    Ny = N0(id,2) + Uy(:,ii);
    plot(Nx,Ny,'Linewidth',1.5,...
        'Color',col(1,ii/size(Ux,2)));
    hold on;
end
axis equal;

function Dist = PneuNet(P,W,H,E,T)
R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,-W/2,E,T,H+H/2);
R3 = dRectangle(P,W-E,W+W/2,T,H+H/2);
C1 = dCircle(P,0,T + 0.5,1);
C2 = dCircle(P,W,T + 0.5,1);
Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end
clr;
%% set signed distance function
W = 20;   % width
H = 40;   % height  
R = 30;   % height  
E = 1;    % edge 
T = 20;   % thickness

sdf = @(x) PneuNet(x,W,H,R);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,W,0,H],...
              'NElem',500);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('MaxIterationMMA',100,...
              'VolumeInfill',0.3,...
              'Penal',4,...
              'FilterRadius',2,...
              'Periodic',[1/2, 0],...
              'Repeat',[],...
              'Nonlinear',false,...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Line',[0,0,0,T]); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Line',[W,W,0,T]); 
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindElements('Location',[W/2,H/2],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[W/2,H/2],.5);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();

function Dist = PneuNet(P,W,H1,H2)
R1 = dRectangle(P,0,W,0,H1);
L1 = dLine(P,W,0,H2,H1);
Dist = dIntersect(R1,L1);
end
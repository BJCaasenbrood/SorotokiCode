clr;
%% generate mesh from sdf
P = 5*kpa;    % pressure
W = 80;         % width
T = 2;          % thickness

Material = Dragonskin10(10);%YeohMaterial('C1',0.0145,'C2',0.00001);

sdf = @(x) dRectangle(x,0,W,0,T);

msh = Mesh(sdf,'BdBox',[0,W,0,T],'Quads',[70 1]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/400);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Location',[0,T/2]),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Location',[W,T/2]),[1,1]);

id = fem.FindEdges('EdgeSelect',[W/2,0],30);
fem = fem.AddConstraint('Pressure',id,[P,0]);

%% select material
fem.Material = Material;

%% solving
fem.solve();

clr;
%% generate mesh from sdf
P = 1*kpa;    % pressure
W = 80;         % width
T = 2;          % thickness

Material = Ecoflex0050(0.1);%YeohMaterial('C1',0.0145,'C2',0.00001);

sdf = @(x) dRectangle(x,0,W,0,T);

msh = Mesh(sdf,'BdBox',[0,W,0,T],'Quads',[70 4]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Location',[0,T/2]),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Location',[W,T/2]),[1,1]);

id = fem.FindEdges('EdgeSelect',[W/2,0],30);
fem = fem.AddConstraint('Pressure',id,[P,0]);

%% select material
fem.Material = Material;

%% solving
fem.solve();

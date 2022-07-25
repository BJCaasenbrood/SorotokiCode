clr;
%% generate mesh from sdf
P = 3*kpa;  % pressure
W = 80;     % width
T = 2;      % thickness

Material = Ecoflex0030(1);

sdf = sRectangle(0,W,0,T);

msh = Mesh(sdf,'BdBox',[0,W,0,T],'Quads',[160 1]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/160,'Linestyle','none','SigmoidFactor',0.5);

%% add constraint
fem = fem.AddConstraint('Spring',fem.FindNodes('SE'),[1,1]*1);
fem = fem.AddConstraint('Spring',fem.FindNodes('SW'),[1,1]*1);
fem = fem.AddConstraint('Spring',fem.FindNodes('NE'),[1,1]*1e3);
fem = fem.AddConstraint('Spring',fem.FindNodes('NW'),[1,1]*1e3);

id = fem.FindEdges('EdgeSelect',[W/2,0],30);
fem = fem.AddConstraint('Pressure',id,P);

%% select material
fem.Material = Material;

%% solving
fem.solve();

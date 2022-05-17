clr;
%% generate mesh from sdf
P = 5*kpa;    % pressure
W = 80;       % width
T = 1.5;      % thickness

Material = Ecoflex0030(150);

sdf = sRectangle(0,W,0,T);

msh = Mesh(sdf,'BdBox',[0,W,0,T],'Quads',[160 1]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/160,'Nonlinear',true,'SigmoidFactor',0.95);

%% add constraint
fem = fem.AddConstraint('Spring',fem.FindNodes('SE'),[1,1]*1e-1);
fem = fem.AddConstraint('Spring',fem.FindNodes('SW'),[1,1]*1e-1);
fem = fem.AddConstraint('Spring',fem.FindNodes('NE'),[1,1]*1e-1);
fem = fem.AddConstraint('Spring',fem.FindNodes('NW'),[1,1]*1e-1);
% fem = fem.AddConstraint('Support',fem.FindNodes('NE'),[1,0]);
% fem = fem.AddConstraint('Support',fem.FindNodes('NW'),[1,0]);
% fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
% fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[0,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Location',[0,T/2]),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Location',[W,T/2]),[1,1]);

id = fem.FindEdges('EdgeSelect',[W/2,0],30);
fem = fem.AddConstraint('Pressure',id,P);

%% select material
fem.Material = Material;

%% solving
fem.solve();

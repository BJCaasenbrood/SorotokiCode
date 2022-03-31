clr;
%% generate mesh from sdf
sdf = sRectangle(0,40,0,2);

msh = Mesh(sdf,'Quads',[50 2]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/120,'SolverPlot',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Right'),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);

%% select material
fem.Material = NeoHookeanMaterial(0.1,0.4);%Ecoflex0030(50);

%% solving
fem.simulate();

%% plotting
%fem.show();
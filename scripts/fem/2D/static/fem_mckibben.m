clr;
%% generate mesh from sdf
P = 1*kpa;  % pressure
W = 15;     % width
H = 80;     % heigth
T = 3;      % thickness
O = 0;      % horizontal offsett

Material = Ecoflex0030(20);

sdf = sRectangle(0,W,0,H) - sRectangle(T-O,W-T-O,5,H-5);
msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',500);
msh = msh.generate();
msh.show;

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/70,'Linestyle','none',...
    'ResidualNorm',1e-1,'Penal',4,'OptimizationProblem','Compliant');

%% add constraint
fem = fem.addSupport(fem.FindNodes('Top'),[1,1]);
fem = fem.addSupport(fem.FindNodes('Bottom'),[1,0]);
fem = fem.addPressure(fem.FindEdges('Hole'),P);
fem = fem.addOutput(fem.FindNodes('Bottom'),[0,1]);

%% select material
fem.Material = Material;

%% solving
fem.solve();

%%
fig(106); clf;
p  = fem.Log.t*P/kpa;
id = fem.FindNodes('Bottom');

dV = fem.Log.Volume - fem.Log.Volume(1);
dX = fem.Log.Out.f;

sorocolor
yyaxis left
plot(p,dV,'LineW',2);

yyaxis right
plot(p,-dX,'LineW',2);

grid on;
box on;
set(gca,'LineW',1.5);
axis square
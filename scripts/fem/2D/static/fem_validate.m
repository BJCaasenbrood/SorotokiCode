clr;
%%
L  = 10;
h  = 1;
b  = 1;
F  = 5;  
E  = 1e6;
Nu = 0;

%%
msh = Mesh(sRectangle(0,L,-h/2,h/2),'Quads',[50,5]);
msh = msh.generate();

msh.show;

%% 
fem = Fem(msh,'TimeStep',1/25,'ResidualNorm',1e-6);
fem.Material = NeoHookeanMaterial(4.125*E/2,Nu);

fem = fem.addSupport('Left',[1,1]);
a = numel(fem.FindNodes('Right'));
fem = fem.addLoad('Right',[0, F / (b * h)]*(1/a));

fem = fem.solve();

%% theory
I = b * h^3 / 12;
w = F * L^3 / (3 * E * I)

wfem = mean(fem.Node(fem.FindNodes('Right'),2))
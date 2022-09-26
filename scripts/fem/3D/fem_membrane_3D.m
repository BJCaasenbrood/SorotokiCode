%% setup 
P0 = 10*kpa; % pressure
W  = 10;     % width in(mm)
T  = 1;    % thickens in (mm)
Ne = 12;     % #top elements

%% building mesh
msh = Mesh(sCube(0,W,0,W,0,T),'Hexahedron',ceil([Ne,Ne,2]));
msh = msh.generate();

%% compute force per elemental face
A0 = W*W/(Ne^2);  % area per element
F  = P0*A0;       % force per node

%% build fem 
fem = Fem(msh,'TimeStep',1/75,'ResidualNorm',1e-4);
fem = fem.addSupport(fem.FindNodes('Left'),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Right'),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Front'),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Back'),[1,1,1]);

fem = fem.addTendon(fem.FindNodes('Bottom'),[0,0,F]);

fem.Material = NeoHookeanMaterial(0.05,0.3);
fem.solve();

%%
figure(101);
fem.replay('axis',[-2,12,-2,12,0,10]);
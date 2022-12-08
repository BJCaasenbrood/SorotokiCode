%% setup 
P0 = 10*kpa; % pressure
W  = 10;     % width in(mm)
T  = 0.5;    % thickens in (mm)
Ne = 15;     % #top elements

%% building mesh
msh = Mesh(sCube(0,W,0,W,0,T),'Hexahedron',ceil([Ne,Ne,2]));
msh = msh.generate();

%% compute force per elemental face
A0 = W*W/(Ne^2);  % area per element
F  = P0*A0;       % force per node

%% build fem 
fem = Fem(msh,'TimeStep',1/75);
fem = fem.addSupport(fem.FindNodes('Left'),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Right'),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Front'),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Back'),[1,1,1]);

fem = fem.addTendon(fem.FindNodes('Bottom'),[0,0,F]);

fem.Material = NeoHookeanMaterial(0.1,0.33);
fem.solve();

%%
figure(101);
fem.replay('axis',[0,10,0,10,0,8]);

id = fem.FindNodes('Location',[W/2,W/2,T]);
for ii = 1:numel(fem.Log.t)
    Y(ii) = fem.Log.Node{ii}(id,3);
end

figure(103)
plot(fem.Log.t*P0/kpa,Y - Y(1),'LineW',2);
xlabel('Pressure (kPa)');
ylabel('Displacement center (mm)');
grid on;
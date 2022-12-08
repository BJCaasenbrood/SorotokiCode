clr;
W = 150;
H = 12;

%% generate mesh
msh = Mesh('Crawler.png','BdBox',[0,W,0,H],...
           'SimplifyTol',0.1,'Hmesh',[1,1,2]);

msh = msh.generate();

%% FEM
fem = Fem(msh,'TimeStep',1/1250,'TimeEnd',6);

fem.Material = NeoHookeanMaterial(0.2,0.4);
fem.Material.Zeta = .75;
fem.Material.Rr   = .3;
fem.Material.Cfr  = .5;

fem = fem.addGravity([0,-9.81e3]);
fem = fem.addContact(SDF(W));

%% make inextensible bottom layer
id = fem.FindElements('Box',[0,150,0,2]);
fem.Density(id) = 5;

%%
fem = fem.addPressure(fem.FindEdges('BoxHole',[0 50 0 15]),...
   @(x) pset(x.Time,1));
fem = fem.addPressure(fem.FindEdges('BoxHole',[50 100 0 15]),...
   @(x) 0.75*pset(x.Time,2));
fem = fem.addPressure(fem.FindEdges('BoxHole',[100 150 0 15]),...
   @(x) pset(x.Time,3));

%% simulate
fem = fem.simulate();

function y = pset(t,id)
w  = 2*pi;
P0 = 15 * kpa;
t  = mod(t,.75);
y  = P0*tsin(w*t - (id-1)*pi/5).*sigmoid(3*w*t);
end

function y = SDF(W)
y = sLine(3*W,-0.1*W,.01,.01);
end
clr;
W = 150;
H = 12;

%% generate mesh
msh = Mesh('Crawler.png','BdBox',[0,W,0,H],...
           'SimplifyTol',0.1,'Hmesh',[1,1,2],...
           'MatlabMeshType','Linear');

msh = msh.generate();
msh.show();

%% FEM
fem = Fem(msh,'TimeStep',1/1250,'TimeEnd',6);

fem.Material = NeoHookeanMaterial(0.2,0.4);
fem.Material.Zeta = .5;
fem.Material.Rr   = .1;
fem.Material.Cfr  = .95;

fem = fem.addGravity([0,-9.81e3]);
fem = fem.addContact(SDF(W));

%%
id = fem.FindElements('Box',[0,150,0,2]);
fem.Density(id) = 5;

%%
fem = fem.addPressure(fem.FindEdges('BoxHole',[0 50 0 15]),...
   @(x) pset(x.Time,1));
fem = fem.addPressure(fem.FindEdges('BoxHole',[50 100 0 15]),...
   @(x) 0.75*pset(x.Time,2));
fem = fem.addPressure(fem.FindEdges('BoxHole',[100 150 0 15]),...
   @(x) pset(x.Time,3));

fem = fem.simulate();

%%
figure(106)
t = linspace(0,6,1e3).';
p = [pset(t,1),pset(t,2),pset(t,3)];

plot(t,p)

%% movie
t = fem.Log.t; close all;
figure(101);

for ii = 1:fps(t,120):numel(t)
    subplot(2,1,1);
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',0*fem.Log.Stress{ii});
    plot([150,150],[0,25],'k','LineW',1.5);
    caxis([-1 1]);
    axis([-1.2*W,1.2*W,-15,5*H]);
    
    subplot(2,1,2);
    P = [pset(t(1:ii).',1).',pset(t(1:ii).',2).',pset(t(1:ii).',3).'];
    plot(t(1:ii).',P/kpa,'LineW',1.5);
    xlim([0, fem.get('TimeEnd')]); xlabel('time (s)','interpreter','latex')
    ylim([0, 5]); ylabel('$p_1$, $p_2$, $p_3$ (kPa)','interpreter','latex')
    
    drawnow();
end

function y = pset(t,id)
w  = 2*pi;
P0 = 15 * kpa;
t = mod(t,.75);
y  = P0*tsin(w*t - (id-1)*pi/5).*sigmoid(3*w*t);
end

function y = SDF(W)
y = sLine(3*W,-0.1*W,.01,.01);
end
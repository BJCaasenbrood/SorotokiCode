clr;
W = 150;
H = 15;

%% generate mesh
msh = Mesh('Crawler.png','BdBox',[0,W,0,H],...
           'SimplifyTol',0.02,'Hmesh',[1,3,5],...
           'MatlabMeshType','quadratic');

msh = msh.generate();
msh.show();

%% FEM
fem = Fem(msh,'TimeStep',1/1250,'TimeEnd',2,...
    'BdBox',[-1.2*W,0.1*W,-15,H],'Linestyle','none');

fem.Material = NeoHookeanMaterial(0.3,0.4);
% fem.Material.Rho = fem.Material.Rho;
fem.Material.Cfr = 5e-6;

fem = fem.addGravity([0,-9.81e3]);
fem = fem.addContact(SDF(W));

%%
fem = fem.addPressure(fem.FindEdges('BoxHole',[0 50 0 15]),...
   @(x) pset(x.Time,1));
fem = fem.addPressure(fem.FindEdges('BoxHole',[50 100 0 15]),...
   @(x) pset(x.Time,2));
fem = fem.addPressure(fem.FindEdges('BoxHole',[100 150 0 15]),...
   @(x) pset(x.Time,3));

fem = fem.simulate();

%%
p = [pset(t,1),pset(t,2),pset(t,3)];


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
w  = 10*pi;
P0 = 5*kpa;
y  = P0*tsin(w*t - (id-1)*pi/2.75).*sigmoid(w*t);
end

function y = SDF(W)
y = sLine(1*W,-0.05*W,-.1,-.10);
end
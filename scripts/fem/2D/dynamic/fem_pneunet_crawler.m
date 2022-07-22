clr;
W = 150;
H = 15;
w  = 5*pi*2;
P0 = 4*kpa;
%% generate mesh
msh = Mesh('Crawler.png','BdBox',[0,W,0,H],...
           'SimplifyTol',0.02,'Hmesh',[1,3,5]);

msh = msh.generate();

msh.show();

%% FEM
fem = Fem(msh,'TimeStep',1/800,'TimeEnd',2.75,...
    'BdBox',[-1.2*W,1.2*W,-15,H],'Linestyle','none');

fem.Material = NeoHookeanMaterial(0.25,0.4);
fem.Material.Rho  = fem.Material.Rho;
fem.Material.Cfr  = 5e-6;

fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Contact',@(x) SDF(x,W),[0,0]);

%%
fem = fem.AddConstraint('Pressure',fem.FindEdges('BoxHole',[0 50 0 15]),...
   @(x) P0*tsin(w*x.Time)*sigmoid(w*x.Time));
fem = fem.AddConstraint('Pressure',fem.FindEdges('BoxHole',[50 100 0 15]),...
   @(x) P0*tsin(w*x.Time - pi/3)*sigmoid(w*x.Time));
fem = fem.AddConstraint('Pressure',fem.FindEdges('BoxHole',[100 150 0 15]),...
   @(x) P0*tsin(w*x.Time - 2*pi/3)*sigmoid(w*x.Time));

fem = fem.simulate();

%% movie
t = fem.Log.t; close all;
figure(101);

for ii = 1:fps(t,300):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    %caxis([0 1e-5]);
    axis([-1.2*W,1.2*W,-15,5*H]);
    drawnow();
end



function D = SDF(x,W)
    %S = sRectangle(-W,2*W,-150,-5);
    S = sLine(3*W,-W,-.1,-.10);
    D = S.eval(x);
end
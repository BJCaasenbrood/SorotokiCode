clr;
%% pressure settings
P0  = 10*kpa;

%% generate mesh
msh = Mesh('PneunetFineDual.png','BdBox',[0,120,-20,20],...
           'SimplifyTol',0.02,'Hmesh',[1,2,3]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/120,'BdBox',[0,120,-80,20],...
              'Linestyle','none','TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);

fem = fem.AddConstraint('Pressure',fem.FindEdges('BoxHole',[0,120,0,20]),...
   @(x) P0*sigmoid(x.Time));

%% assign material
fem.Material = NeoHookeanMaterial(0.65,0.3);
fem.Material.Zeta = 0.75;

%% solve
fem.solve(); 

%% extract dominated modes
NNode = 200;
Modal = [0,2,0,1,0,0];
X     = linspace(0,1,NNode)';

% shape function reconstruction
shp = Shapes(fem,Modal,'NNode',NNode,...
    'L0',120,'FilterRadius',[1,5]);

shp.g0 = SE3(eye(3),[0;0;1e-1]);
shp = shp.reference([0,0],[119,0]);
shp = shp.rebuild();
shp = shp.reconstruct();

shp.show();

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;
h = [];

shp.Kp = 1e-2;
shp.Kd = 1e6;

for ii = 1:fps(t,50):numel(t)
    delete(h);
    
    Nds = fem.Log.Node{ii};
    fem.set('Node',Nds);
    fem.show('Field',fem.Log.Stress{ii});
    axis([0 130 -110 110]);
    
    Xi = shp.recoverStrain(fem,ii);
    q  = shp.estimateJointSpace(Xi,Nds);
    p  = shp.FK(q);
    
    h = fplot(p(:,[1,3]),'LineW',3,'Color',magenta);
    
    background();
    drawnow();
end

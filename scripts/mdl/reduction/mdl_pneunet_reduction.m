clr;
%% pressure settings
P0  = 10*kpa;

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.02,'Hmesh',[1,6,12]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/35,'BdBox',[0,120,-80,20],...
              'Linestyle','none','TimeEnd',2.5);
        
fem.Material = NeoHookeanMaterial(0.65,0.3);
fem.Material.Zeta = 1.75;

%% solve quasi-static          
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem.solve();

Uqs = fem.Log.U(end,:).';

%% adding pressure loads
fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time)*(0.5 + 0.5*sin(x.Time)));

%% solve dynamics
fem.set('Utmp',Uqs);
fem.simulate(); 

%% extract dominated modes
NNode = 200;
Modal = [0,3,0,1,0,0];
X     = linspace(0,1,NNode)';

% shape function reconstruction
shp = Shapes(fem,Modal,'NNode',NNode,...
    'L0',120,'FilterRadius',[15,15]);

shp.g0 = SE3(eye(3),[0;0;1e-1]);
shp = shp.reference([0,1e-1],[119,1e-1]);
shp = shp.reconstruct();

shp.show();

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;
h = [];

shp.Kp = 1;
shp.Kd = 1e8;

Q = [];

for ii = 1:numel(t)
    delete(h);
    
    Nds = fem.Log.Node{ii};
    fem.set('Node',Nds);
    fem.show('Field',fem.Log.Stress{ii});
    axis([-50 120 -100 30]);
    
    Xi = shp.recoverStrain(fem,ii);
    q  = shp.estimateJointSpace(Xi,Nds);
    p  = shp.FK(q);
    
    h = fplot(p(:,[1,3]),'LineW',3,'Color',magenta);
    
    Q = [Q; q.'];
    
    background();
    drawnow();
end

%%
figure(102);
plot(t,Q);
%plot(Q,gradient(Q));


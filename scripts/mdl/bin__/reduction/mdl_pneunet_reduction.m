clr;
%% pressure settings
w   = 0.1;
P0  = 20*kpa;

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.02,'Hmesh',[1,2,2.5],...
           'MatlabMeshType','quadratic');

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/15,'BdBox',[0,120,-80,20],...
              'Linestyle','none','TimeEnd',2.5);
        
fem.Material = NeoHookeanMaterial(0.65,0.4);%Dragonskin10(25);
fem.Material.Zeta = 1.75;

%% solve quasi-static          
fem = fem.addSupport(fem.FindNodes('Box',[0,0,0,10]),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.addGravity([9.81e3,0]);
fem.solve();

Uqs = fem.Log.U(end,:).';

%% adding pressure loads
fem = fem.addGravity([9.81e3,0]);
fem = fem.addPressure(fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time));

%% solve dynamics
fem.set('Utmp',Uqs,'TimeStep',1/60,'ResidualNorm',1e-2);
fem.simulate(); 

%% extract dominated modes
NNode = 200;
Modal = [0,3,0,1,0,0];
X     = linspace(0,1,NNode)';

% shape function reconstruction
shp = Shapes(fem,Modal,'NNode',NNode,...
    'L0',120,'FilterRadius',[10,10]);

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

shp.Kp = 5;
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


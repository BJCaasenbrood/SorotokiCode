clr;
%% pressure settings
P0  = 5*kpa;

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.02,'Hmesh',[1,2,3],...
           'MatlabMeshType','quadratic');

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/120,'BdBox',[0,120,-80,20],...
              'Linestyle','none','TimeEnd',2);

%% add boundary constraint
fem = fem.addSupport(fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.addGravity([0,-9.81e3]);

fem = fem.addPressure(fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time));

%% assign material
fem.Material = NeoHookeanMaterial(0.65,0.3);
fem.Material.Zeta = 0.75;

%% solve
fem.simulate(); 

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,120):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    axis([-50 120 -100 30]);
    
    background();
    drawnow();
end

function Dist = SDF(x,R)
    sdf  = sCircle(30,-50 + R,R);
    Dist = sdf.eval(x);
end
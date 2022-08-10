clr;
%% pressure settings
P0  = 15*kpa;

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.03,'Hmesh',[1,2,3],...
           'MatlabMeshType','linear');

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/300,'BdBox',[-20,120,-80,20],...
              'Linestyle','none','TimeEnd',2);

%% add boundary constraint
fem = fem.addSupport(fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.addContact(@(x) SDF(x,10,70,-30),[0,0]);
fem = fem.addGravity([0,-9.81e3]);

fem = fem.addPressure(fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time));
%% assign material
fem.Material = NeoHookeanMaterial(0.35,0.3);
fem.Material.Zeta = 1.45;

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
    axis([-5 130 -40 60]);
    
    background();
    drawnow();
end

function Dist = SDF(x,R,x0,y0)
    sdf  = sCircle(x0,y0 + R,R);
    Dist = sdf.eval(x);
end
clr;
%% pressure settings
P0  = 25*kpa;

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.03,'Hmesh',[1,2,3]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/50,'BdBox',[-20,120,-80,20],...
              'Linestyle','none','TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
%fem = fem.AddConstraint('Contact',@(x) SDF(x,15), [0, 0]);

fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time));

%% assign material
fem.Material = NeoHookeanMaterial(0.65,0.3);
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
    axis([-50 120 -100 30]);
    
    background();
    drawnow();
end

function Dist = SDF(x,R)
    sdf  = sCircle(10,-40 + R,R);
    Dist = sdf.eval(x);
end
clr;
%% pressure settings
P0  = 20*kpa;

%% generate mesh
msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.02,'Hmesh',[1,2,4]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/250,'BdBox',[0,120,-80,20],...
              'Linestyle','none','TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
%fem = fem.AddConstraint('Contact',@(x) SDF(x,22),[0,0]);

fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time));

%% add output nodes
% id  = fem.FindNodes('Bottom');
% fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = Dragonskin10();
fem.Material.Zeta = 1.4;

%% solve
fem.simulate(); 

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,320):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    axis([-30 120 -100 30]);
    
    background(gitpage);
    drawnow();
end

function Dist = SDF(x,R)
sdf  = sCircle(30,-50 + R,R);
Dist = sdf.eval(x);
end
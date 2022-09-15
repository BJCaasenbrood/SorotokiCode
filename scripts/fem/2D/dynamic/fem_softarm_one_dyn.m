clr;
%% generate mesh
Simp  = 0.08;
GrowH = 1;

MinH  = 6;
MaxH  = 8;

msh = Mesh('SWR_robot_1.png',...
    'BdBox',[0,250,0,290],...
    'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh.show();
msh.show('Holes');

%% pressure sets
P0    = 3*kpa;
pset0 = [23,21,18,15,12,9,6,3];

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/50,...
              'BdBox',[0,108,0,325],...
              'Linestyle','none',...
              'MovieAxis',[-25 120 -60 130],...
              'Movie',0,'TimeEnd',4);

%% add boundary constraint
fem = fem.addSupport(fem.FindNodes('Bottom'),[1,1]);
fem = fem.addPressure(fem.FindEdges('Hole',...
    pset0), @(x) P0*sigmoid(x.Time));

%% assign material
fem.Material = NeoHookeanMaterial(2,0.25);  
fem.Material.Zeta = 0.35;

%% solve
fem.solve(); 

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,20):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    axis([0,350,0,400]);
    background();
    drawnow();

end


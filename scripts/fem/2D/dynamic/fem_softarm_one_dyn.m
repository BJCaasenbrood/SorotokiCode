clr;
%% generate mesh
Simp  = 0.08;
GrowH = 1;

MinH  = 5;
MaxH  = 8;

msh = Mesh('SWR_robot_1.png','BdBox',[0,250,0,290],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh.show();
msh.show('Holes');

%% pressure sets
P0    = 3*kpa;
pset0 = [23,21,18,15,12,9,6,3];

%pset0 = [25,21,18,15,12,9,7,4];
%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/50,'BdBox',[0,108,0,325],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',4);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('Hole',pset0), @(x) P0*sigmoid(x.Time));

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
    subplot(1,3,[1,2]);
        fem.set('Node',fem.Log.Node{ii});
        fem.show('Field',fem.Log.Stress{ii});
        axis([-300,300,-325,325]);      
        
    subplot(1,3,3);
        cla;
        K = fem.Log.Kin(1:ii);
        plot(t(1:ii),K,'LineW',2); hold on;
        %axis([0 1 -8e-3 7e-3]);
       
        U0 = fem.Log.Vg(end) + fem.Log.Psi(end) - fem.Log.Vf(end);
        U = fem.Log.Psi(1:ii) + fem.Log.Vg(1:ii) - fem.Log.Vf(1:ii) - U0;
        plot(t(1:ii),U,'LineW',2); 
        axis([0 t(end) -1 3]);     

        plot(t(1:ii),K+U,'LineW',2); 
    
    drawnow();
end


clr;
%% generate mesh
Simp  = 0.025;
GrowH = 2;

MinH  = 1;
MaxH  = 3;

%msh = Mesh('SWR_robot_gripper.png','BdBox',...
msh = Mesh('SWR_robot.png','BdBox',...
    [0,108,0,325],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();
msh.show();
msh.show('Holes');

%% pressure sets
P0 = -1.5*kpa;
P1 = 1.5*kpa;

pset0 = [0 100 250 400];
pset1 = [0 40 0 120];
pset2 = [60 100 0 120];
pset3 = [0 40 120 240];
pset4 = [60 100 120 240];

%pset0 = [25,21,18,15,12,9,7,4];
%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/125,'BdBox',[0,120,0,325],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',3);

%% add boundary constraint
fem = fem.addSupport(fem.FindNodes('Box',[0 150 0 2]),[1,1]);
fem = fem.addGravity([0,-1e3]);

fem = fem.addPressure(fem.FindEdges('BoxHole',pset1), @(x) 2*P1*sigmoid(x.Time));
fem = fem.addPressure(fem.FindEdges('BoxHole',pset2), @(x) 0.75*P0*sigmoid(x.Time));
fem = fem.addPressure(fem.FindEdges('BoxHole',pset3), @(x) 0.75*P0*sigmoid(x.Time));
fem = fem.addPressure(fem.FindEdges('BoxHole',pset4), @(x) 2*P1*sigmoid(x.Time));

%% assign material
fem.Material = NeoHookeanMaterial(2,0.33);  
fem.Material.Zeta = 0.35;

%% solve
fem.simulate(); 

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,30):numel(t)
    %subplot(1,3,[1,2]);
        fem.set('Node',fem.Log.Node{ii});
        fem.show('Field',fem.Log.Stress{ii});
        axis([-300,300,-325,325]);      
        
%     subplot(1,3,3);
%         cla;
%         K = fem.Log.Kin(1:ii);
%         plot(t(1:ii),K,'LineW',2); hold on;
%         %axis([0 1 -8e-3 7e-3]);
%        
%         U0 = fem.Log.Vg(end) + fem.Log.Psi(end) - fem.Log.Vf(end);
%         U = fem.Log.Psi(1:ii) + fem.Log.Vg(1:ii) - fem.Log.Vf(1:ii) - U0;
%         plot(t(1:ii),U,'LineW',2); 
%         axis([0 t(end) -1 3]);     
% 
%         plot(t(1:ii),K+U,'LineW',2); 
    
    drawnow();
end

%% Signed distance environment
function D = SDF(x)
    %S = sLine(130,130,300,0);
    S = sCircle(53,262,16);
    D = S.eval(x);
end
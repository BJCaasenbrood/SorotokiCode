clr;
%% pressure settings
P0  = 1*kpa;

%% generate mesh
Simp  = 0.06;
GrowH = 1;

MinH  = 20;
MaxH  = 30;

msh = Mesh('SWR_robot.png','BdBox',[0,624,0,1604],...
    'SimplifyTol',Simp,'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

figure(101); clf;
msh.show();

%% re-orient the mesh
H1 = [29,32,35,38,41,44,46,48];
H2 = [28,31,34,37,40,43,47,50];
H3 = [4,7,10,13,16,19,22,25];
H4 = [3,6,9,12,15,18,21,24];

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/30,'BdBox',[0,120,-80,20],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
%fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);

fem = fem.AddConstraint('Pressure',fem.FindEdges('Hole',H1),...
    @(t) -0.2*P0*sigmoid(t));

fem = fem.AddConstraint('Pressure',fem.FindEdges('Hole',H2),...
    @(t) P0*sigmoid(t));

fem = fem.AddConstraint('Pressure',fem.FindEdges('Hole',H3),...
    @(t) P0*sigmoid(t));

fem = fem.AddConstraint('Pressure',fem.FindEdges('Hole',H4),...
    @(t) -0.3*P0*sigmoid(t));
% 
% %% add output nodes
% id  = fem.FindNodes('Bottom');
% fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = Dragonskin30(10);

%% solve
% out = input('* Simulate again? (y/n) ','s');
% if strcmp(out,'y')
fem.solve(); 
%     close all;
%     save('job.mat');
% else
%     load('job.mat'); 
% end

%% shape function reconstruction
clearvars -except fem;
shp = Shapes(fem,[0,1,0,0,0,0],'NModal',[2],'NNode',400);
shp = shp.reference([320,30],[320,1550]);
shp = shp.reconstruct();

%%
figure(102); clf;

subplot(3,2,[1,2]);
X   = shp.Sigma;
Phi = shp.POD{1};

plot(shp.Sigma,Phi/abs(max(Phi)),'LineW',3,'Col',col(1)); hold on;
Phi = shp.POD{2};

plot(shp.Sigma,Phi/abs(max(Phi)),'LineW',3,'Col',col(2));
axis on; set(gca,'LineW',1.5);
grid on; box on;
axis([0 1500 -2 2]);

k = 1;
for ii = [0,-.5,-1,-1.5,-2,-2.5]
subplot(3,2,[3,5]); hold on;
[p, g] = shp.string([ii,0]); 
plot(p(:,1),-p(:,3),'LineW',3,'Col',col(k));
%axis([-5 135 -120 10]);
axis equal;
axis on; set(gca,'LineW',1.5);
grid on; box on;
k = k +1;
end
legend([])

k = 1;
for ii = [0,0.05,0.1,0.15,0.2,0.25]
subplot(3,2,[4,6]); hold on;
[p, g] = shp.string([0,ii]); 
plot(p(:,1),-p(:,3),'LineW',3,'Col',col(k));
%axis([-5 125 -120 10]);
axis equal;
axis on; set(gca,'LineW',1.5);
grid on; box on;
k = k +1;
end

%plot(p(:,1),-p(:,3),'LineW',2);
%plot(p(:,1),-p(:,3)*0,'--k','LineW',1);
%axis([-30 150 -150 20]); hold on;
%axis equal;

%% post-processing
figure(106); %cla;
[p, g] = shp.string([-1,0]); hold on;
plot(p(:,1),-p(:,3),'LineW',3);
% t = fem.Log.t; close all;
% figure(105);
% 
% for ii = 1:fps(t,10):numel(t)
%     N = fem.Log.Node{ii};
%     %subplot(1,2,1);
%     fem.set('Node',N);
%     fem.show('Field',fem.Log.Stress{ii});
%     axis([-500 1000 0 2000]);
%     background(gitpage);
%     drawnow();
% end
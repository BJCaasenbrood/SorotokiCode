clr;
%% pressure settings
P0  = 30*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;

MinH  = 2;%0.8586;
MaxH  = 3;%1.5;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

figure(101);
subplot(2,1,1); imshow('Pneunet.png');
subplot(2,1,2); msh.show();
%msh.NElem

%% re-orient the mesh
%msh = Blender(msh,'Rotate',{'x',90}); msh = msh.show(); pause(2);

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/125,'BdBox',[0,120,-80,20],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
    @(t) P0*sigmoid(t));
   % @(t) P0*(0.4 + 0.6*sin(2*t^2))*min(2*t,1));

%% add output nodes
id  = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = NeoHookeanMaterial(1,0.25);%Dragonskin30(30);

%% solve
out = input('* Simulate again? (y/n) ','s');
if strcmp(out,'y')
    clf;
    tic;
    [~,tmpNode] = fem.simulate(); 
    toc;
    close all;
    save('job.mat');
else
    load('job.mat'); 
end

%% shape function reconstruction
clearvars -except fem shp;
shp = Shapes(fem,[0,1,0,1,0,0],'NModal',[2,1],'NNode',40);
shp = shp.reconstruct();
save('job.mat');
%%
figure(102); clf;

% subplot(3,2,[1,2]);
% X   = shp.Sigma;
% Phi = shp.POD{1};
% plot(X/max(X),Phi/abs(max(Phi)),'LineW',3,'Col',col(1)); hold on;
% Phi = shp.POD{2};
% plot(X/max(X),Phi/abs(max(Phi)),'LineW',3,'Col',col(2));
% axis on; set(gca,'LineW',1.5);
% grid on; box on;
% axis([0 1 -2 2]);

k = 1;
for ii = [0,-.5,-1,-1.5,-2,-2.5]
subplot(3,2,[3,5]); hold on;
[p, g] = shp.string([ii,0]); 
plot(p(:,1),-p(:,3),'LineW',3,'Col',colfade(k,5));
axis equal;
axis on; set(gca,'LineW',1.5);
grid on; box on;
axis([-30 145 -100 10]);
k = k +1;
end

k = 1;
for ii = [0,0.05,0.1,0.15,0.2,0.25]
subplot(3,2,[4,6]); hold on;
[p, g] = shp.string([0,ii]); 
plot(p(:,1),-p(:,3),'LineW',3,'Col',colfade(k,5));

axis equal;
axis on; set(gca,'LineW',1.5);
grid on; box on;
axis([-30 145 -100 10]);
k = k +1;
end

%% figure POD
figure(107); clf;
X   = shp.Sigma;
subplot(3,2,[1,3])
for ii = 1:3
    Phi = shp.POD{ii};
    plot(X/max(X),Phi/abs(max(Phi)),'LineW',3,'Col',col(ii)); hold on;
end
axis on; set(gca,'LineW',1.5);
grid on; box on;
axis([0 1 -2 2]);

subplot(3,2,[2,4])
for ii = 1:3
    Phi = shp.PODQ{ii};
    plot(X/max(X),Phi/abs(max(Phi)),'LineW',3,'Col',col(ii)); hold on;
end
axis on; set(gca,'LineW',1.5);
grid on; box on;
axis([0 1 -2 2]);

subplot(3,2,[5]);
plot(shp.PODEnergy{1}/max(shp.PODEnergy{1}),'-o','LineW',3);
axis on; set(gca,'LineW',1.5); grid on;
axis([0,10,0,1]);

subplot(3,2,[6]);
plot(shp.PODEnergy{2}/max(shp.PODEnergy{2}),'-o','LineW',3);
axis on; set(gca,'LineW',1.5); grid on;
axis([0,10,0,1]);

%% energy balance
t  = fem.Log.t;
K  = fem.Log.Kin;
Ue = fem.Log.Psi;
Ug = fem.Log.Vg;
Uf = fem.Log.Vf;

figure(108); cla;
plot(t,K+Ue+Ug-Uf,'LineW',2,'Col',col(1)); 
%plot(t,Ue); hold on;
%plot(t,Ug); hold on;
%% post-processing
% t = fem.Log.t; close all;
% figure(105);
% 
% for ii = 1:fps(t,30):numel(t)
%     N = fem.Log.Node{ii};
%     %subplot(1,2,1);
%     fem.set('Node',N);
%     fem.show('Field',fem.Log.Stress{ii});
%     axis([-60 130 -120 30]);
%     background(gitpage);
%     drawnow();
% end

function C = colfade(ii,X)
C = col(1,min((ii)/(X),1));
end
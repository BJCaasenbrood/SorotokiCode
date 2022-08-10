clr;
%% pressure settings
T  = 4;
DW = 0.5*pi;
N  = 15;
t  = linspace(0,T,1e4);

figure(99); clf;
for ii = 1:N
    W = DW*ii;
    subplot(1,2,1);
    plot(t,yd(t,W),'LineW',1.5); hold on;
    
    subplot(1,2,2);
    plot(t,gradient(yd(t,W)),'LineW',1.5); hold on;
end

%pause();

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.02,'Hmesh',[1,2,3],...
           'MatlabMeshType','quadratic');

msh = msh.generate();
msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/35,'BdBox',[0,120,-80,20],...
              'Linestyle','none','TimeEnd',T,'ShowProcess',1);
        
fem.Material = NeoHookeanMaterial(0.65,0.4);
fem.Material.Zeta = 1.75;

%% solve quasi-static          
fem = fem.addSupport(fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.addGravity();
fem.solve();

Uqs = fem.Log.U(end,:).';
fem = fem.set('Utmp',Uqs);

k = 1;
for ii = (N):-1:1
    W = DW*ii;
    fem = fem.set('Utmp',Uqs);
    fem = fem.set('TimeStep',1/(100));
    
    disp(['Excitation: w = ',num2str(W)]);
    try
        POD(:,:,k) = evalFrequency(fem,W);
        disp(['success']);
    catch me
        disp(['fail']);
    end
    k = k + 1;
end

%%
figure(102);
clf;
color = evolution(6);
for kk = 1:5
    subplot(1,5,kk);
for ii = 1:5
    plot(POD(:,kk,ii),'Color',color(ii,:),'LineW',2);
    hold on;
    axis square;
end
end

%%
function POD = evalFrequency(fem,w)
% adding pressure loads
fem = fem.addPressure(fem.FindEdges('AllHole'),...
    @(x) yd(x.Time,w));

% solve dynamics
fem.simulate();

% extract dominated modes
NNode = 200;
Modal = [0,3,0,2,0,0];

% shape function reconstruction
shp = Shapes(fem,Modal,'NNode',NNode,...
    'L0',120,'FilterRadius',[15,15]);

shp.g0 = SE3(eye(3),[0;0;1e-1]);
shp = shp.reference([0,1],[119,1]);
shp = shp.reconstruct();

POD = shp.get('PODR');
end

function y = yd(t,w)
 P0 = 5*kpa;
 y  = P0*sigmoid(0.5*t,0.5).*(0.5 + 0.5*sin(w*t));
end

%shp.show();

% %% movie
% close all;
% fem.set('Linestyle','none');
% figure(101);
% 
% t = fem.Log.t;
% h = [];
% 
% shp.Kp = 5;
% shp.Kd = 1e8;
% 
% Q = [];
% 
% for ii = 1:numel(t)
%     delete(h);
%     
%     Nds = fem.Log.Node{ii};
%     fem.set('Node',Nds);
%     fem.show('Field',fem.Log.Stress{ii});
%     axis([-50 120 -100 30]);
%     
%     Xi = shp.recoverStrain(fem,ii);
%     q  = shp.estimateJointSpace(Xi,Nds);
%     p  = shp.FK(q);
%     
%     h = fplot(p(:,[1,3]),'LineW',3,'Color',magenta);
%     
%     Q = [Q; q.'];
%     
%     background();
%     drawnow();
% end

%%
%figure(102);
%plot(t,Q);
%plot(Q,gradient(Q));


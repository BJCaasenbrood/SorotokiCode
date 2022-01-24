%clr; cdsoro; beep off;
%clc; clear;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',8,'NDisc',1);
mdl = setupSoftRobot(mdl,0.01,0.001,0.001);
%mdl = setupSoftRobot(mdl,1.0,0.1,0.1);
mdl = mdl.set('Controller',1);
mdl = mdl.set('TimeStep', 1/33.33333);
mdl = mdl.set('Tdomain', 10); 

%% generate and solve dynamic model
mdl = mdl.generate();
mdl = mdl.csolve(); 

%% show results                                               
figure(102)
t  = mdl.get('t');

% subplot(3,4,[1 2 5 6]);
subplot(2,1,1);
% plot(t,mdl.q_,'k-','linewidth',0.5); hold on; grid on;
plot(t,mdl.q,'linewidth',2);  grid on;
%     
% subplot(3,4,[3 4 7 8]);
% plot(t,mdl.H,'linewidth',1.5);hold on; grid on;
% legend('$\mathcal{T}$','$\mathcal{V}_e$','$\mathcal{V}_g$','$\mathcal{H}$',...
%     'interpreter','latex','FontSize',12);
% 
% subplot(3,4,9:12);
% plot(t,mdl.ge(:,5:end),'linewidth',1.5); hold on; grid on;
% plot(t,mdl.gd(:,5:end),'k--','linewidth',1); 
% legend('$x$','$y$','$z$',...
%     'interpreter','latex','FontSize',12);
% 
% figure(102)
% subplot(2,1,1);
% plot(t,mdl.q,'linewidth',1.5);  grid on;
% ylabel('Modal coefficients $q(t)$','interpreter','latex','fontsize',12)
% 
subplot(2,1,2);
hold on;
plot(t,mdl.ge(:,5:end),'linewidth',2); hold on; grid on;
plot(t,mdl.gd(:,5:end),'k--','linewidth',1.5); 
% box on;
% legend('$x$','$y$','$z$',...
%     'interpreter','latex','FontSize',12);
% ylabel('Modal coefficients $q(t)$','interpreter','latex','fontsize',12)
%plot(t,mdl.H(:,4),'linewidth',1.5);hold on; grid on;
%plot(t(2:end),diff(mdl.H(:,4))/mean(diff(t)),'linewidth',1.5);hold on; grid on;
%plot(t,mdl.H(:,4),'linewidth',1.5);hold on; grid on;
%plot(t,mdl.detae(:,2),'linewidth',1.5);hold on; grid on;
% plot(t,mdl.H(:,4) + mdl.detae(:,1),'linewidth',1.5);hold on; grid on;

% subplot(2,1,2);
% hold on;
% plot(t,mdl.H(:,1) ,'linewidth',1.5);
% %plot(t,mdl.H(:,4)- mdl.detae(:,1),'linewidth',1.5);hold on; grid on;
% box on;
% legend('$H_s(q,p)$',...
%    'interpreter','latex','FontSize',12);

%plot(t(2:end),diff(mdl.H(:,4) + mdl.detae(:,1)),'linewidth',1.5);hold on; grid on;


grid on;
%error(1);
disp(' - Press enter to play simulation - ');
%% generate rig
[rig, sph] = setupRig(mdl);

%text(0.055,0.00,-0.005,'\textbf{$g_d$}','interpreter','latex','fontsize',16);

h = [];
for ii = 1:fps(t,15):length(mdl.q)
    rig = rig.computeFK(mdl.q(ii,:));
    rig = rig.update();
   
    sph.reset();
    sph = Blender(sph,'SE3x',mdl.gd(ii,:));
    sph = Blender(sph,'SE3',rig.g0);
    sph.update();

    setupFigure(ii);
    view(30,30);
    
    delete(h);
    %h = shadowplot(5);
    
    background();
    
%     if mod(ii,15) <= 1
%         pause;
%         t(ii)
%     end
%     if ii == 1, gif('srm3_octarm.gif','frame',gcf,'nodither');
%         pause; framepause(5);
%     else, gif;
%     end
end

%framepause(15);

%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,Kp,Kd,Lam)
L0 = 0.12;

mdl = mdl.set('Sdomain',   L0);
mdl = mdl.set('SpaceStep', 24);
mdl = mdl.set('Density',   1200);
mdl = mdl.set('Radius',    1e-2);
mdl = mdl.set('Gravity',   [0,0,-9.81]);
mdl = mdl.set('E',         25);
mdl = mdl.set('Mu',        0.2);
mdl = mdl.set('Gain',      [Kp,Kd,0]);
mdl = mdl.set('Spring',    [0.01e-7,1]);
mdl = mdl.set('Lambda',    [Lam,0]);%[Kp/Lam,Kp*30*0]);

mdl = mdl.set('ActuationSpace',-1);
mdl = mdl.set('Movie',1);

mdl = mdl.set('Point',...
    [1,0,0,0,0.04,0.00,0.005]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl = Gmodel('Arm.stl');

gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
    'SSSPower',0.1,'SSSRadius',0.15,'SSS',true);

rig = Rig(@(x) mdl.string(x),'Domain',mdl.Sdomain);
rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,.999);

rig = rig.texture(1,mateplastic);
rig.g0 = [rot2quat(roty(pi/2)),0,0,0];

sph = Gmodel('Sphere.stl');
sph = Blender(sph,'Scale',3e-3);
sph.Texture = 1.25*diffuse(0.95);
sph = sph.fix();
sph = sph.render();
sph.update();

rig = rig.render();
end

% setup figure
function setupFigure(ii)
axis([0 0.12 -0.06 0.06 -0.06 0.04]);
box on;
axis off;
if (ii == 1)
    view(30,30); 
    drawnow;
else
    drawnow;
end
background(metropolis);
end

function framepause(k)

i = 1;
while i < k
   gif
   i = i+1;
end
end

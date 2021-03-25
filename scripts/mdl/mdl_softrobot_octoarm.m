clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',8,'NDisc',1);
mdl = setupSoftRobot(mdl,10e-5,70,500);
mdl = mdl.set('TimeStep', 1/50);
mdl = mdl.set('Tdomain',  10); 

%% generate and solve dynamic model
mdl = mdl.generate();
mdl.q0(1) = -7;
mdl.q0(2) = 15;
mdl = mdl.csolve(); 

%% show results                                               
figure(102)
t  = mdl.get('t');

subplot(3,4,[1 2 5 6]);
plot(t,mdl.q_,'k-','linewidth',0.5); hold on; grid on;
plot(t,mdl.q,'linewidth',1.5);  grid on;

subplot(3,4,[3 4 7 8]);
plot(t,mdl.H,'linewidth',1.5);hold on; grid on;
legend('$\mathcal{T}$','$\mathcal{V}_e$','$\mathcal{V}_g$','$\mathcal{H}$',...
    'interpreter','latex','FontSize',12);

subplot(3,4,9:12);
plot(t,mdl.ge(:,5:end),'linewidth',1.5); hold on; grid on;
plot(t,mdl.gd(:,5:end),'k--','linewidth',1); 
legend('$x$','$y$','$z$',...
    'interpreter','latex','FontSize',12);

disp(' - Press enter to play simulation - ');
%pause;
%% generate rig
[rig, sph] = setupRig(mdl);

%text(0.06,0,0.025,'\textbf{$g_d$}','interpreter','latex','fontsize',16);

for ii = 1:fps(t,7):length(mdl.q)
    rig = rig.compute(ii);
    rig = rig.update();
   
    sph.reset();
    sph = Blender(sph,'SE3',mdl.gd(ii,:));
    sph = Blender(sph,'SE3',rig.g0);
    sph.update();

    setupFigure(ii);
    view(0,0);
end

%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,Kp,Kd,Lam)
mdl = mdl.set('Controller',1);
L0 = 0.12;

mdl = mdl.set('Sdomain',   L0);
mdl = mdl.set('SpaceStep', 150);
mdl = mdl.set('Density',   800);
mdl = mdl.set('Radius',    0.01);
mdl = mdl.set('Gravity',   [0,0,-9.81]);
mdl = mdl.set('E',         25);
mdl = mdl.set('Mu',        0.1);
mdl = mdl.set('Gain',      [Kp,Kd,Kp*5e3]);
mdl = mdl.set('Spring',    [1e-3,1]);
mdl = mdl.set('Lambda',    [Kp/Lam,Kp*30]);

mdl = mdl.set('ActuationSpace',-1);
mdl = mdl.set('Movie',1);

mdl = mdl.set('Point',...
    [1,0,0,0,0.06,0.00,-0.01]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl = Gmodel('Arm.stl');

gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
    'SSSPower',0.1,'SSSRadius',0.15,'SSS',true);

rig = Rig(mdl);
rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,.999);

rig = rig.texture(1,mateplastic);
rig.g0 = [rot2quat(roty(-pi/2)).',0,0,0];

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
axis([0 0.1 -0.06 0.06 -0.08 0.05]);
box on;
axis off;
if (ii == 1)
    view(30,30); 
    drawnow;
else
    drawnow;
end
end


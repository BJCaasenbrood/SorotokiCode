clr;
%% assign free DOF
mdl = Model([0,1,1,1,0,0],'NModal',1,'NDisc',1);
mdl = setupSoftRobot(mdl,50e-5,20,500);

mdl = mdl.set('TimeStep',1e-3);
mdl = mdl.set('Tdomain', 3); 

%% generate and solve dynamic model
mdl = mdl.generate(); mdl.q0(1) = 5;
mdl = mdl.csolve(); 

%% show results                                               
figure(102)
t  = mdl.get('t');

subplot(3,4,[1 2 5 6]);
plot(t,mdl.etae,'linewidth',1.5); hold on; grid on;

subplot(3,4,[3 4 7 8]);
plot(t,mdl.H,'linewidth',1.5);hold on; grid on;
legend('$\mathcal{T}$','$\mathcal{V}_e$','$\mathcal{V}_g$','$\mathcal{H}$',...
    'interpreter','latex','FontSize',12);

subplot(3,4,9:12);
plot(t,mdl.ge(:,5:end),'linewidth',1.5); hold on; grid on;
plot(t,mdl.gd(:,5:end),'k--','linewidth',1); 
legend('$x$','$y$','$z$',...
    'interpreter','latex','FontSize',12);

%% generate rig
[rig, sph] = setupRig(mdl);

for ii = 1:fps(t,30):length(mdl.q)
    rig = rig.compute(ii);
    rig = rig.update();
      
    sph.reset();
    sph = Blender(sph,'SE3',mdl.gd(ii,:));
    sph = Blender(sph,'SE3',rig.g0);
    sph.update();

    setupFigure(ii);
    title(['T = ',num2str(t(ii),3)]);
end

%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,Kp,Kd,Lam)
mdl = mdl.set('Controller',0);
L0 = 0.12;

mdl = mdl.set('Sdomain',   L0);
mdl = mdl.set('SpaceStep', 25);
mdl = mdl.set('Density',   1800);
mdl = mdl.set('Radius',    0.01);
mdl = mdl.set('Gravity',   [0,0,0]);
mdl = mdl.set('E',         15e2);
mdl = mdl.set('Mu',        0.05);
mdl = mdl.set('Gain',      [Kp,Kd]);
mdl = mdl.set('Spring',    [1e-10,1]);
mdl = mdl.set('Lambda',    Kp/Lam);

mdl = mdl.set('ActuationSpace',1);

mdl = mdl.set('Point',...
    [1,0,0,0,0.065,0.00,-0.025]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl = Gmodel('Pneulink.stl');
gmdl.set('TextureStretch',.85);

rig = Rig(mdl);
rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,.975);

rig.g0 = [rot2quat(roty(-pi)).',0,0,0];

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
axis([-0.05 0.05 -0.05 0.05 -0.1 0]);
box on;
axis off;
if (ii == 1)
    view(30,30); 
    drawnow;
else
    drawnow;
end
end


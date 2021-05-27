clr; cdsoro;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',6,'NDisc',2);
mdl = setupSoftRobot(mdl,1e3,3,20);
mdl = mdl.set('Controller',1);
mdl = mdl.set('TimeStep',1/150);
mdl = mdl.set('Tdomain',20); 

%% generate trajectory
mdl = mdl.generate(); 

%% generate and solve dynamic model
mdl = mdl.csolve(); 

%% show results
figure(102);

t  = mdl.get('t');

subplot(3,4,[1 2 5 6]);
plot(t,mdl.q,'linewidth',1.5); hold on; grid on;

subplot(3,4,[3 4 7 8]);
plot(t,mdl.H,'linewidth',1.5);hold on; grid on;
legend('$\mathcal{H} + \mathcal{H}_c$',...
    'interpreter','latex','FontSize',12);

subplot(3,4,9:12);
plot(t,mdl.ge(:,5:end),'linewidth',1.5); hold on; grid on;
plot(t,mdl.gd(:,5:end),'k--','linewidth',1); 
legend('$x$','$y$','$z$',...
    'interpreter','latex','FontSize',12);
%error('bla');
%% generate rig
[rig, sph] = setupRig(mdl);

for ii = 1:fps(t,7):length(mdl.q)
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
mdl = mdl.set('Controller',1);

L0 = 0.36;

mdl = mdl.set('Sdomain',   L0);
mdl = mdl.set('SpaceStep', 50);
mdl = mdl.set('Density',   150);
mdl = mdl.set('Radius',    0.05);
mdl = mdl.set('Gravity',   [0,0,0]);
mdl = mdl.set('E',         250);
mdl = mdl.set('Mu',        0.05);
mdl = mdl.set('Gain',      [Kp,Kd,0]);
mdl = mdl.set('Spring',    [0,1]);
mdl = mdl.set('Lambda',    [Kp/Lam,0]);
mdl = mdl.set('ActuationSpace',0);

mdl = mdl.set('Point',...
    [0,0,0,0,0.275,0.1,0.05]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl1 = Gmodel('Pneulink.stl');
gmdl2 = Gmodel('Pneulink.stl');
gmdl3 = Gmodel('SoftGripperRedux.stl');
gmdl3 = Blender(gmdl3,'Scale',0.475);

assignin('base','gmdl1',gmdl1);
assignin('base','gmdl2',gmdl2);
assignin('base','gmdl3',gmdl3);

rig = Rig(mdl);
rig = rig.add(gmdl1,gmdl2,gmdl3);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,0.475);
rig = rig.parent(2,0,0.475);
rig = rig.parent(2,1,0.95);
rig = rig.parent(3,0,0.95);

rig = rig.texture(1:3,base*1.15);
rig.g0 = [rot2quat(roty(pi)).',0,0,0];

sph = Gmodel('Sphere.stl');
sph = Blender(sph,'Scale',15e-3);
sph.Texture = diffuse(0.95);
sph = sph.fix();
sph = sph.render();

axis([-0.07 0.07 -0.07 0.07 -0.15 0]);
view(-30,10);
rig = rig.render();
end

% setup figure
function setupFigure(ii)
axis(3*[-0.06 0.06 -0.06 0.06 -0.12 0]);
box on;
axis off;
if (ii == 1)
    view(-30,10);
    drawnow;
else
    drawnow;
end

end



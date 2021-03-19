clr;
%% assign free DOF
mdl = Model([0,1,0,1,0,0],'NModal',4,'NDisc',1);
mdl = setupSoftRobot(mdl,0.5,0,1e2);
mdl = mdl.set('Tdomain',15); 
mdl = mdl.set('Movie',0); 

%% generate and solve dynamic model
mdl = mdl.generate();
mdl.q0_(1) = 15;
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

%% generate rig
% figure(106);
% clf;
% subplot(2,1,1); 
% plot(t,mdl.etae(:,1:3),'linewidth',1.5); 
% %plot(t,mdl.detae,'linewidth',1.5); 
% %axis([0 10 -.02 .02])
% % 
% subplot(2,1,2); hold on;
% plot(t,mdl.etae(:,4:6),'linewidth',1.5); 
% axis([0 10 -.2 .2]);

if 1
[rig, sph] = setupRig(mdl);

text(0.05,0,-0.06,'\textbf{$g_d$}','interpreter','latex','fontsize',16);

for ii = 1:fps(t,8):length(mdl.q)
    
    rig = rig.compute(ii);
    rig = rig.update();
      
    sph.reset();
    sph = Blender(sph,'SE3',mdl.Point);
    sph = Blender(sph,'SE3',rig.g0);
    sph.update();
   
    setupFigure(ii);
    title(['T = ',num2str(t(ii),3)]);
    
    if ii == 1, pause(0.5); end
    mdl = mdl.updateFrame();
end
end

%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,Kp,Kd,Lam)
mdl = mdl.set('Controller',1);

L0 = 0.063;
mdl = mdl.set('TimeStep',   1/50);
mdl = mdl.set('Sdomain',    L0);
mdl = mdl.set('SpaceStep',  10);
mdl = mdl.set('Density',    50);
mdl = mdl.set('Radius',     0.03);
mdl = mdl.set('Gravity',    [0,-9.81,0]);
mdl = mdl.set('E',          250);
mdl = mdl.set('Mu',         0.05);
mdl = mdl.set('Gain',       [Kp,Kd,Kp]);
mdl = mdl.set('Lambda',     [Kp/Lam,Kp*1e-2]);


mdl = mdl.set('ActuationSpace',1);

[px,py] = SpherePosition(-15,L0*1.2);

mdl = mdl.set('Point',...
    [0,0,0,0,py,0,px]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl1 = Gmodel('SoftActuatorPlanarRedux.stl');
gmdl1.set('TextureStretch',.85);
assignin('base','gmdl1',gmdl1);

rig = Rig(mdl);
rig = rig.add(gmdl1);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,0.99);

rig = rig.texture(1,base);
rig.g0 = [rot2quat(roty(pi)).',0,0,0];

sph = Gmodel('Sphere.stl');
sph = Blender(sph,'Scale',5e-3);
sph.Texture = diffuse(0.93);
sph = sph.fix();

sph = sph.render();
rig = rig.render();

axis([-0.05 0.05 -0.03 0.03 -0.1 0]);
view(30,10); 

assignin('base','sph',sph);
end

% compute possible desired trajectory
function [px,py] = SpherePosition(k,l)
    a = 1/(sign(k)*abs(k) + 1e-3);
    px = a*(1 - cos(l*k));
    py = a*sin(l*k);
end

% setup figure
function setupFigure(ii)
axis([-0.05 0.05 -0.03 0.03 -0.1 0]);
box on;
axis off;
if (ii == 1)
    view(30,10); 
    drawnow;
else
    drawnow;
end
end


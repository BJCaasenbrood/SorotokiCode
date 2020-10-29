clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',6,'NDisc',1);
mdl = setupSoftRobot(mdl,1e-5,1.0);

%% generate and solve dynamic model
mdl = mdl.generate();
mdl = mdl.csolve(); 

%% show results
figure(102)
t = mdl.get('t');
ge = mdl.ge;
xd = mdl.get('xd');

subplot(3,4,[1 2 5 6]);
shade(t,mdl.q,'linewidth',1.5); hold on;

subplot(3,4,[3 4 7 8]);
plot(t,ge(:,5:end),'linewidth',1.5); hold on;
plot(t,xd(:,5:end),'k--','linewidth',1.0); hold on;
% 
subplot(3,4,9:12);
plot(t,mdl.get('tau'),'linewidth',1.0); 
%% generate rig
[rig, sph] = setupRig(mdl);

for ii = 1:fps(t,6):length(mdl.q)
    rig = rig.compute(ii);
    rig = rig.update();
      
    sph.reset();
    sph = Blender(sph,'SE3',xd(ii,:));
    sph = Blender(sph,'SE3',rig.g0);
    sph.update();

    setupFigure(ii);
    title(['T = ',num2str(t(ii),3)]);
end

%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,K,X)
mdl = mdl.set('Controller', 1);
L0 = 0.12;

mdl = mdl.set('Tdomain',    20); 
mdl = mdl.set('TimeStep',   1/12);
mdl = mdl.set('Sdomain',    L0);
mdl = mdl.set('SpaceStep',  25);
mdl = mdl.set('Density',    50);
mdl = mdl.set('Radius',     0.01);
mdl = mdl.set('Gravity',    [0,0,-9.81]);
mdl = mdl.set('E',          400);
mdl = mdl.set('Mu',         0.2);
mdl = mdl.set('Gain',       [K,0]);
mdl = mdl.set('Lambda',     K/X);

mdl = mdl.set('ActuationSpace',-1);

mdl = mdl.set('Point',...
    [1,0,0,0,0.06,0.02,-0.02]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl = Gmodel('Arm.stl');

gmdl.set('Emission', [0.9 0.8 0.8],...
    'SSSPower',0.1,'SSSRadius',0.15,'SSS',true);

rig = Rig(mdl);
rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig = rig.texture(1,mateplastic);
rig.g0 = [rot2quat(roty(-pi/2)).',0,0,0];

sph = Gmodel('Sphere.stl');
sph = Blender(sph,'Scale',4e-3);
sph.Texture = 1.25*diffuse(0.6);
sph = sph.fix();
sph = sph.render();

rig = rig.render();
end

% setup figure
function setupFigure(ii)
axis([0 0.1 -0.05 0.05 -0.05 0.05]);
box on;
axis on;
if (ii == 1)
    view(60,30); 
    drawnow;
else
    drawnow;
end
end


clr;
%% assign free DOF
mdl = Model([0,1,0,0,0,0],'NModal',5,'NDisc',1);
mdl = setupSoftRobot(mdl,1e-4,0.75,25);

%% generate and solve dynamic model
mdl = mdl.generate();
mdl.q0(1) = 1e-3;
mdl = mdl.csolve(); 

%% show results
figure(105)
t = mdl.get('t');
q = mdl.q;
ge = mdl.ge;
u = mdl.get('tau');
xd = mdl.get('xd');

subplot(3,4,[1 2 5 6]);
plot(t,(q),'linewidth',1.0); hold on;

subplot(3,4,[3 4 7 8]);
plot(t,ge(:,5:end),'linewidth',1.0); hold on;
plot(t,xd(:,5:end),'k--','linewidth',1.0); hold on;

subplot(3,4,9:12);
plot(t,u,'linewidth',1.0); 

%% generate rig
rig = setupRig(mdl);

sph = Gmodel('Sphere.stl');
sph = Blender(sph,'Scale',5e-3);
sph.Texture = prusa;
sph = sph.fix();
sph = sph.render();

FPS = 2;

for ii = 1:FPS:length(mdl.q)
    rig = rig.compute(ii);
    rig = rig.update();
    
    sph.reset();
    sph = Blender(sph,'SE3',xd(ii,:));
    sph.update();

    axis([-0.05 0.05 -0.05 0.05 .009 0.1]);
end

%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,Kp,Kd,R)
L0 = 0.15;

mdl = mdl.set('Controller',1);

mdl = mdl.set('Tdomain',    10); 
mdl = mdl.set('TimeStep',   1/12);
mdl = mdl.set('Sdomain',    L0);
mdl = mdl.set('SpaceStep',  50);
mdl = mdl.set('Density',    2);
mdl = mdl.set('Radius',     0.02);
mdl = mdl.set('Gravity',    9.81);
mdl = mdl.set('E',          15);
mdl = mdl.set('Mu',         0.2);
mdl = mdl.set('Gain',       [Kp,0]);
mdl = mdl.set('Lambda',     Kp/Kd);

[px,py] = SpherePosition(R,L0);

mdl = mdl.set('Point',...
    [1,0,0,0,py,0,px]);
end

% setup rig
function rig = setupRig(mdl)
gmdl = Gmodel('Pneunet.stl');

assignin('base','gmdl',gmdl);

rig = Rig(mdl);
rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);
rig = rig.texture(1,diffuse(100));
rig = rig.render();

rig.g0 = [rot2quat(roty(-pi/2)).',0,0,0];

figure(101);
gmdl.ground([-0.04 0.02 -0.03 0.03 -0.00 0.1]);

axis([-0.05 0.05 -0.05 0.05 .009 0.1]);
end

% compute possible desired trajectory
function [px,py] = SpherePosition(k,l)
    a = 1/(sign(k)*abs(k) + 1e-3);
    px = a*(1 - cos(l*k));
    py = a*sin(l*k);
end



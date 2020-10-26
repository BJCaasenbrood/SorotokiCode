clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',5,'NDisc',1);
mdl = setupSoftRobot(mdl,1e-7,0.25);

%% generate and solve dynamic model
mdl = mdl.generate();
mdl = mdl.csolve(); 

%% show results
figure(102)
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
%% generate rig
rig = setupRig(mdl);

sph = Gmodel('Sphere.stl');
sph = Blender(sph,'Scale',2e-3);
sph.Texture = prusa;
sph = sph.fix();
sph = sph.render();

FPS = 3;

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
function mdl = setupSoftRobot(mdl,K,X)
mdl = mdl.set('Controller',1);

mdl = mdl.set('Tdomain',    25); 
mdl = mdl.set('TimeStep',   1/12);
mdl = mdl.set('Sdomain',    0.12);
mdl = mdl.set('SpaceStep',  50);
mdl = mdl.set('Density',    500);
mdl = mdl.set('Radius',     0.02);
mdl = mdl.set('Gravity',    9.81);
mdl = mdl.set('E',          50);
mdl = mdl.set('Mu',         0.1);
mdl = mdl.set('Gain',       [K,0]);
mdl = mdl.set('Lambda',     K/X);

mdl = mdl.set('Point',...
    [1,0,0,0,0.05,0.025,0.025]);
end

% setup rig
function rig = setupRig(mdl)
gmdl = Gmodel('Arm.stl');
gmdl.Alpha = 1.0;

mus1 = Gmodel('Cylinder.stl');
mus1 = Blender(mus1,'Scale',{'axi',0.015});
mus1 = Blender(mus1,'Rotate',{'y',2.5});
mus1 = Blender(mus1,'Translate',{'x',0.05});
mus1 = mus1.fix();

mus2 = Gmodel('Cylinder.stl');
mus2 = Blender(mus2,'Scale',{'axi',0.015});
mus2 = Blender(mus2,'Rotate',{'y',2.5});
mus2 = Blender(mus2,'Translate',{'x',0.05});
mus2 = Blender(mus2,'Rotate',{'z',180});
mus2 = mus2.fix();

mus3 = Gmodel('Cylinder.stl');
mus3 = Blender(mus3,'Scale',{'axi',0.015});
mus3 = Blender(mus3,'Rotate',{'y',2.5});
mus3 = Blender(mus3,'Translate',{'x',0.05});
mus3 = Blender(mus3,'Rotate',{'z',90});
mus3 = mus3.fix();

mus4 = Gmodel('Cylinder.stl');
mus4 = Blender(mus4,'Scale',{'axi',0.015});
mus4 = Blender(mus4,'Rotate',{'y',2.5});
mus4 = Blender(mus4,'Translate',{'x',0.05});
mus4 = Blender(mus4,'Rotate',{'z',-90});
mus4 = mus4.fix();

assignin('base','gmdl',gmdl);
assignin('base','mus1',mus1);
assignin('base','mus2',mus2);
assignin('base','mus4',mus4);

rig = Rig(mdl);
rig = rig.add(gmdl,mus1,mus2,mus3,mus4);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);
rig = rig.parent(2,0,0.1);
rig = rig.parent(2,1,0.95);
rig = rig.parent(3,0,0.1);
rig = rig.parent(3,1,0.95);
rig = rig.parent(4,0,0.1);
rig = rig.parent(4,1,0.9);
rig = rig.parent(5,0,0.1);
rig = rig.parent(5,1,0.9);
rig = rig.texture(1,redwax);
rig = rig.texture(2:3,redgloss);
rig = rig.texture(4:5,metal);
rig = rig.hide(2,3,4,5);
rig = rig.render();
gmdl.ground([-0.015 0.015 -0.015, 0.015 0]);

axis([-0.05 0.05 -0.05 0.05 .009 0.1]);
end



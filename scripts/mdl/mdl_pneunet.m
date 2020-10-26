clr;
%% assign free DOF
mdl = Model([0,1,0,0,0,0],'NModal',5,'NDisc',1);
mdl = setupSoftRobot(mdl,1e-7,0.25);

%% generate and solve dynamic model
mdl = mdl.generate();
mdl.q0(1) = -20;
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
FPS = 3;

for ii = 1:FPS:length(mdl.q)
    rig = rig.compute(ii);
    rig = rig.update();
    

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
mdl = mdl.set('Density',    20);
mdl = mdl.set('Radius',     0.02);
mdl = mdl.set('Gravity',    -9.81);
mdl = mdl.set('E',          50);
mdl = mdl.set('Mu',         0.1);
mdl = mdl.set('Gain',       [K,0]);
mdl = mdl.set('Lambda',     K/X);

mdl = mdl.set('Point',...
    [1,0,0,0,0.09,0.025,0]);
end

% setup rig
function rig = setupRig(mdl)
gmdl = Gmodel('Pneunet.stl');

assignin('base','gmdl',gmdl);

rig = Rig(mdl);
rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);
rig = rig.texture(1,base);
rig = rig.render();

gmdl.ground([-0.05 0.05 -0.05 0.05 .009 0.1]);

axis([-0.05 0.05 -0.05 0.05 .009 0.1]);
end



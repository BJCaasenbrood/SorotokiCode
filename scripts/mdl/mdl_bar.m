clr;
%% assign free DOF
mdl = Model([1,1,0,0,0,0],'NModal',5,'NDisc',1);
mdl = setupSoftRobot(mdl,1e-3);

%% generate and solve dynamic model
mdl = mdl.generate();
mdl.q0 = 1e-3*rand(mdl.NDof*mdl.NModal,1);
mdl = mdl.csolve(); 

%% show results
figure(102)
t = mdl.get('t');
ge = mdl.ge;
xd = mdl.get('xd');
P = mdl.get('Phi');

subplot(3,4,[1 2 5 6]);
plot(t,mdl.q,'linewidth',1.5); hold on; grid on;

subplot(3,4,[3 4 7 8]);
plot(t,ge(:,5:end),'linewidth',1.5); hold on;
plot(t,xd(:,5:end),'k--','linewidth',1.0); hold on;
grid on;

%% generate rig
rig = setupRig(mdl);

for ii = 1:fps(t,6):length(mdl.q)
    rig = rig.compute(ii);
    rig = rig.update();
      
    setupFigure(ii);
    title(['T = ',num2str(t(ii),3)]);
end



%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,K)
mdl = mdl.set('Controller', -1);
L0 = 1.0;
W = 10e-3;
H = 40e-3;

mdl = mdl.set('Tdomain',    40); 
mdl = mdl.set('TimeStep',   1e-3);
mdl = mdl.set('Sdomain',    L0);
mdl = mdl.set('SpaceStep',  12);
mdl = mdl.set('Radius',     0.01);
mdl = mdl.set('E',          50);
mdl = mdl.set('Mu',         4);
mdl = mdl.set('Lambda',     K);


mdl = mdl.set('Jxx', (W*H^2)*(1/12));
mdl = mdl.set('Jyy', (H*W^2)*(1/12));
mdl = mdl.set('Jzz', (1/12)*(W^2 + H^2));

mdl = mdl.set('ActuationSpace',-1);

mdl = mdl.set('Point',...
    [rot2quat(rotx(-0.2*pi)*roty(-pi))',0.5,0.0,-0.1]);
end

% setup rig
function rig = setupRig(mdl)
gmdl = Gmodel('SlenderRod.stl');

gmdl = Blender(gmdl,'Scale',{'y',4});
assignin('base','gmdl',gmdl);

rig = Rig(mdl);
rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);
rig = rig.texture(1,greyresin);

rig.g0 = [rot2quat(roty(-pi/2)).',0,0,0];

rig = rig.render();
end

% setup figure
function setupFigure(ii)
axis(10*[0 0.1 -0.03 0.03 -0.03 0.03]);
box on;
axis on;
if (ii == 1)
    view(60,30); 
    drawnow;
else
    drawnow;
end
end


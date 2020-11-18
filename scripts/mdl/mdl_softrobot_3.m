clr;
%% assign free DOF
mdl = Model([0,1,0,1,0,0],'NModal',2,'NDisc',1);
mdl = setupSoftRobot(mdl,1,1);

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

subplot(3,4,[1 2 5 6]); overwrite_colors;
plot(t,(q),'linewidth',1.5); hold on;

subplot(3,4,[3 4 7 8]);
plot(t,ge(:,5:end),'linewidth',1.5); hold on;
plot(t,xd(:,5:end),'k--','linewidth',1.0); hold on;
% 
subplot(3,4,9:12);
plot(t,u,'linewidth',1.0); 
%% generate rig
[rig, sph] = setupRig(mdl);

for ii = 1:fps(t,12):length(mdl.q)
    
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

L0 = 0.063;

mdl = mdl.set('Tdomain',    15); 
mdl = mdl.set('TimeStep',   1/30);
mdl = mdl.set('Sdomain',    L0);
mdl = mdl.set('SpaceStep',  50);
mdl = mdl.set('Density',    10);
mdl = mdl.set('Radius',     0.03);
mdl = mdl.set('Gravity',    [-9.81,0,0]);
mdl = mdl.set('E',          50);
mdl = mdl.set('Mu',         0.3);
mdl = mdl.set('Gain',       [K,0.1]);
mdl = mdl.set('Lambda',     K/X);

mdl = mdl.set('ActuationSpace',1);

[px,py] = SpherePosition(15,L0*1.12);

mdl = mdl.set('Point',...
    [1,0,0,0,py,0,-px]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl1 = Gmodel('SoftActuatorPlanarRedux.stl');
gmdl1.set('TextureStretch',.85);
assignin('base','gmdl1',gmdl1);

rig = Rig(mdl);
rig = rig.add(gmdl1);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,0.95);

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
axis on;
if (ii == 1)
    view(30,10); 
    drawnow;
else
    drawnow;
end
end


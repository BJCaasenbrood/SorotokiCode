clr;
%% assign free DOF
mdl = Model([0,0,1,0,0,0],'NModal',5,'NDisc',2);
mdl = setupSoftRobot(mdl,1e-3,0.05);

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
plot(t,(q),'linewidth',1.5); hold on;

subplot(3,4,[3 4 7 8]);
plot(t,ge(:,5:end),'linewidth',1.5); hold on;
plot(t,xd(:,5:end),'k--','linewidth',1.0); hold on;
% 
subplot(3,4,9:12);
plot(t,u,'linewidth',1.0); 

% figure(103); 
% for ii = 1:fps(t,5):length(mdl.q)
%     [u,X] = recoverField(mdl,ii);
%     figure(105);
%     hold on;
%     subplot(2,1,1);
%     plot(X,u(:,1),'Color',...
%         lerp([1,1,1],col(1),ii/200),'linewidth',1.5); 
%     axis([0 0.36 -25 25]);
%     
%     hold on;
%     subplot(2,1,2);
%         plot(X,u(:,2),'Color',...
%         lerp([1,1,1],col(1),ii/200),'linewidth',1.5); 
%   
%     axis([0 0.36 -25 25]);
% end

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

L0 = 0.36;

mdl = mdl.set('Tdomain',    15); 
mdl = mdl.set('TimeStep',   1/12);
mdl = mdl.set('Sdomain',    L0);
mdl = mdl.set('SpaceStep',  30);
mdl = mdl.set('Density',    10);
mdl = mdl.set('Radius',     0.03);
mdl = mdl.set('Gravity',    [-9.81,0,0]);
mdl = mdl.set('E',          250);
mdl = mdl.set('Mu',         0.83);
mdl = mdl.set('Gain',       [K,0.1]);
mdl = mdl.set('Lambda',     K/X);

mdl = mdl.set('ActuationSpace',1);

mdl = mdl.set('Point',...
    [1,0,0,0,0.15,0.05,-0.15]);
end

% setup rig
function [rig, sph] = setupRig(mdl)
gmdl1 = Gmodel('SoftActuatorPlanarRedux.stl');
gmdl2 = Gmodel('SoftActuatorPlanarRedux.stl');
gmdl3 = Gmodel('SoftActuatorPlanarRedux.stl');
%gmdl4 = Gmodel('SoftGripperRedux.stl');

%gmdl4 = Blender(gmdl4,'Scale',0.3);

assignin('base','gmdl1',gmdl1);
assignin('base','gmdl2',gmdl2);
assignin('base','gmdl3',gmdl3);
%assignin('base','gmdl4',gmdl4);

rig = Rig(mdl);
rig = rig.add(gmdl1,gmdl2,gmdl3);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,0.3);
rig = rig.parent(2,0,0.3);
rig = rig.parent(2,1,0.6);
rig = rig.parent(3,0,0.6);
rig = rig.parent(3,1,0.9);
% rig = rig.parent(4,1,0.9);

rig = rig.texture(1:3,base*1.15);
rig.g0 = [rot2quat(roty(pi)).',0,0,0];

sph = Gmodel('Sphere.stl');
sph = Blender(sph,'Scale',15e-3);
sph.Texture = diffuse(0.95);
sph = sph.fix();
sph = sph.render();

axis([-0.07 0.07 -0.07 0.07 -0.15 0]);
view(60,30);
rig = rig.render();
end

% setup figure
function setupFigure(ii)
axis(3*[-0.06 0.06 -0.06 0.06 -0.12 0]);
box on;
axis on
if (ii == 1)
    view(60,30); 
    drawnow;
else
    drawnow;
end

end



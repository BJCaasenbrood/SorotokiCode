clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',6,'NDisc',2);
mdl = setupSoftRobot(mdl);

%% generate and solve dynamic model
mdl = mdl.generate();
mdl = mdl.csolve(); 

%% show sim
figure(102)
t = mdl.get('t');
q = mdl.q;
u = mdl.tau;
ge = mdl.ge;
dq = mdl.dq;
xd = mdl.get('xd');

subplot(3,4,[1 2 5 6]);
plot(t,(q),'linewidth',1.0); hold on;
% 
subplot(3,4,[3 4 7 8]);
plot(t,ge(:,5:end),'linewidth',1.0); hold on;
plot(t,xd(:,5:end),'k--','linewidth',1.0); hold on;
% subplot(3,4,9:12);
% plot(t,u,'linewidth',1.0); 

mdl.showModel();
 
 %% setup model
function mdl = setupSoftRobot(mdl)
mdl = mdl.set('MovieAxis',0.11*[-0.75 0.75 -0.75...
    0.75 -1.75 .1]*0.85,'Movie',0);

mdl = mdl.set('Texture',grey);
mdl = mdl.set('Controller',1);

mdl = mdl.set('Tdomain',    75); 
mdl = mdl.set('TimeStep',   1/12);
mdl = mdl.set('Sdomain',    0.12);
mdl = mdl.set('SpaceStep',  60);
mdl = mdl.set('Density',    12);
mdl = mdl.set('Radius',     0.02);
mdl = mdl.set('Gravity',    9.81);
mdl = mdl.set('E',          15);
mdl = mdl.set('Mu',         0.35);
mdl = mdl.set('Gain',       [1e-4,0.0]);
mdl = mdl.set('Lambda',     5e-4);

mdl = mdl.set('Point',...
    [1,0,0,0,0.1,0,0.04,0]);
 end



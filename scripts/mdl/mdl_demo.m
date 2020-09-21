clr;
%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',4,'NDisc',2);
mdl = mdl.set('MovieAxis',[-0.75 0.75 -0.75 0.75 -1.75 .1]*0.85,'Movie',0);
mdl = mdl.set('Texture',base);
%mdl = mdl.set('Controller',0);

mdl = mdl.set('Point',[1,0,0,0,1,0,0]);

%% generate dynamic model
mdl = mdl.generate();
%mdl.q0(1) = 5;
%% simulate soft robot
mdl = mdl.csolve(); 
% 
%% show sim
figure(102)
t = mdl.get('t');
q = mdl.q;
u = mdl.tau;
ge = mdl.ge;
dq = mdl.dq;
xd = mdl.get('xd');

subplot(3,4,[1 2 5 6]);
plot(t,(q),'linewidth',1.0);
subplot(3,4,[3 4 7 8]);

plot(t,ge(:,5:end),'linewidth',1.0); hold on;
plot(t,xd(:,5:end),'k--','linewidth',1.0); hold on;
subplot(3,4,9:12);
plot(t,u,'linewidth',1.0); 

mdl.showModel();

%% show strain field
figure;

P = mdl.get('Phi');
s = linspace(0,1,101);
Q = mdl.q(end,:);
v = [];

for ii = 1:length(s)
   A = P(s(ii));
   v = vappend(v,(A*Q.').');
end

for jj = 1:size(v,2)    
   subplot(size(v,2),1,jj);
   shade(s,v(:,jj),'Color',col(jj))
end


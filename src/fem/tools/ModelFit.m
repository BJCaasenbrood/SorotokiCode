function Q = ModelFit(P,Anchors,L0,Dofs,varargin)
Npts = 20;
Bpts = 2;

BdBox = boxhull(P);
tol = 0.3*sqrt((BdBox(2)-BdBox(1))*...
    (BdBox(4)-BdBox(3))/size(P,1));

P(:,1) = P(:,1) - mean(P(P(:,2) <= tol,1));
P(:,2) = P(:,2) - mean(P(P(:,2) <= tol,2));

mdl = Model(Dofs,'NModal',6);
mdl = mdl.set('Sdomain', L0);
mdl = mdl.set('SpaceStep',Npts);
mdl = mdl.set('NDisc',1);

for ii = 1:2:length(varargin)
    mdl.set(varargin{ii},varargin{ii+1});
end

mdl = mdl.generate();

C2 = P(2,:);

Q = linspacen([0,0],C2,Bpts+2);
x0 = Q(1:2,2:(2+Bpts-1)); x0 = x0(:);
obj = @(x) ObjectiveBezier(P,[0,0],C2,x);

options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'FiniteDifferenceStepSize',1e-3,'Display','off');

[x, ~, ~, ~] = fminunc(obj,x0,options);

Q = [[0;0],vec2waypoints(x),C2'];
t = linspace(0,1,50);
[C,~,K] = BezierEval(Q,t);

obj = @(x) ObjectiveCosserat(mdl,C',x);
x0 = mdl.q0;

options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'FiniteDifferenceStepSize',1e-4,'Display','off',...
    'FunctionTolerance',1e-2);

[x, ~, ~, ~] = fminunc(obj,x0,options);

figure(103); clf;
plot(C(1,:),C(2,:),'r-','linewidth',2); hold on;
plot(Q(1,:),Q(2,:),'k--o','linewidth',1); hold on;
plot(P(:,1),P(:,2),'b.');

mdl.q0 = x;
g = mdl.string(0,Npts);   
plot(g(:,7),g(:,5),'g-x'); 
axis equal;
end

function f = ObjectiveBezier(P,C1,C2,x)

Q = [C1',vec2waypoints(x),C2'];

t = linspace(0,1,50);
C = BezierEval(Q,t);
[~,D] = distance2curve(C',P,'linear');

f = D.'*D;
end

function f = ObjectiveCosserat(Model,C,x)
Model.q0 = x;
g = Model.string(0,length(C));
G = [g(:,7),g(:,5)];

%f1 = (1/Model.get('E'))*double(x.'*inv(Model.Kee)*x);
f = sqrt(sum((C(:,1)-G(:,1)).^2 + (C(:,2)-G(:,2)).^2));


%[~,D] = distance2curve(G,C,'linear');
%f = D.'*D;
end

function W = vec2waypoints(x)
W = reshape(x,[2,length(x)/2]);
end

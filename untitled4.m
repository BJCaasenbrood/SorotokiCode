clr;
%% extract data
mm2m = 1e-3;
A = readmatrix('data.xlsx');

t = A(:,1);
A(:,1:2:end) = [];

Ux = A(:,1:15);
Uy = A(:,16:end);

%% create spatial data
Sx = linspace(0,70,15);
Sy = linspace(0,70,15)*0;
Px = repmat(Sx,[length(Ux),1]) + Ux;
Py = repmat(Sy,[length(Ux),1]) + Uy;

for ii = 1:3:length(Ux)
    px = mm2m*Px(ii,:);
    py = mm2m*Py(ii,:);

    hold on;
    plot(-py',px','o-','Color',col(1),'linewidth',2);
end

axis([-0.02 0.07 0 0.08]);

%% soft robot model
mdl = Model([0,1,0,1,0,0],'NModal',4,'NDisc',1);
mdl = setupSoftRobot(mdl,50e-5,20,500);
mdl = mdl.generate();

X = [];
x0 = mdl.q0;

figure(108)
for ii = 2:4:length(Ux)
    % optimization
    px = mm2m*Px(ii,:);
    py = mm2m*Py(ii,:);
    P = [py',px'];
    
    obj = @(x) Objective(mdl,P,x);
    
    options = optimoptions('fminunc','Algorithm',...
        'quasi-newton','FiniteDifferenceStepSize',1e-4);
    [x, fval, exitflag, output] = fminunc(obj,x0,options);
    
    %cla;
    % plotting
    mdl.q0 = x;
    g = mdl.string(0,15);
    hold on;
    plot(-py',px','o-','Color',col(1),'linewidth',2);
    plot(-g(:,7),g(:,5),'-o','Color',col(2),'linewidth',1);
    drawnow();
    
    pause(0.5);
    
    X = [X,x];
    x0 = x;
end

%%
s = linspace(0,0.07,50);
P = get(mdl,'Phi');

figure
Y1 = [];
Y2 = [];

for ii = 1:size(X,2)
    Q = X(:,ii);
    
    y = [];
    for jj = 1:length(s)
        xi_ = P(s(jj))*Q;
        y = [y,xi_];
    end
    
    Y1 = [Y1,y(1,:).'];
    Y2 = [Y2,y(2,:).'];
end

figure(103)
cla;
subplot(2,1,1);
for ii = 3:size(Y1,2)
    shade(s/0.07,Y1(:,ii),'Color',col(1),'linewidth',1,'FillAlpha',0.2);
    hold on;
end

subplot(2,1,2);
for ii = 2:size(Y2,2)
    shade(s/0.07,Y2(:,ii),'Color',col(1),'linewidth',1,'FillAlpha',0.2);
    hold on;
end

%% generate POD basis
W1 = Y1*Y1.';
W2 = Y2*Y2.';

[u1,S1,~] = svd(W1);
[u2,S2,~] = svd(W2);

figure(105);
for ii = 1:5
    subplot(2,1,1);
    plot(s/max(s),u1(:,ii),'Color',col(ii),'linewidth',2);
    hold on;
    
    subplot(2,1,2);
    plot(s/max(s),u2(:,ii),'Color',col(ii),'linewidth',2);
    hold on;
end



function f = Objective(Model,P,x)
Model.q0 = x;
G = Model.string(0,15);

f = sqrt(sum((P(:,1)-G(:,7)).^2 + (P(:,2)-G(:,5)).^2));
end

function mdl = setupSoftRobot(mdl,Kp,Kd,Lam)
mdl = mdl.set('Controller',0);
L0 = 0.07;

mdl = mdl.set('TimeStep', 1/5);
mdl = mdl.set('Tdomain',  15); 
mdl = mdl.set('Sdomain',   L0);
mdl = mdl.set('SpaceStep', 15);
mdl = mdl.set('Density',   800);
mdl = mdl.set('Radius',    0.01);
mdl = mdl.set('Gravity',   [0,0,-9.81]);
mdl = mdl.set('E',         25);
mdl = mdl.set('Mu',        0.1);
mdl = mdl.set('Gain',      [Kp,Kd]);
mdl = mdl.set('Spring',    [1e-10,1]);
mdl = mdl.set('Lambda',    Kp/Lam);

mdl = mdl.set('ActuationSpace',-1);

mdl = mdl.set('Point',...
    [1,0,0,0,0.065,0.00,-0.025]);
end





clr; beep off;
%% settings
L = 1;   % manipulator length
M = 3;   % number of modes
N = 301; % grid on SR

%% build basis
Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,0,0,0,0]);          % generate basis

shp.Xi0 = @(x) [0;0;0;1;0;0];
shp = shp.rebuild();

%%
q = zeros(M,1);
q(3) = 7;

%% render curve
p = shp.FK(q);
[~,Et] = shp.tangentPoint(q);

subplot(2,1,1);
cplot(p(2:end,1),p(2:end,3),Et,turbo,'LineW',5);
axis equal;

subplot(2,1,2);
plot(shp.Sigma,Et,'LineW',1.5);

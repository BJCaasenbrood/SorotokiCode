clear; close all; clc;
addpath('build');
addpath('build/log');
addpath(genpath('model'));

%% assign free DOF
mdl = Model([0,1,1,0,0,0],'NModal',5,'NDisc',1);
mdl = setupSoftRobot(mdl,1,.5);

%% generate and solve dynamic model
mdl = mdl.generate();
mdl = mdl.csolve();

figure(103);
subplot(1,2,1);
hold on;
plot(mdl.t,mdl.q,'linewidth',1.5);
subplot(1,2,2);
hold on;
plot(mdl.t,mdl.H,'linewidth',1.5);

figure(102);
for ii = 1:5:length(mdl.t)
    [g, X] = string(mdl,ii);
    figure(102);
    plot3(g(:,5),g(:,6),g(:,7),'b');
    axis equal;
    view(30,30);
    axis(0.2*[-1 1.5 -1 1 -1 1]);
    pause(1/24);
end

%% BACK-END FUNCTIONS
% setup model
function mdl = setupSoftRobot(mdl,K,X)
mdl = mdl.set('Controller', 0);
L0 = 0.25;

mdl = mdl.set('Tdomain',   15); 
mdl = mdl.set('TimeStep',  1/25);
mdl = mdl.set('Sdomain',   L0);
mdl = mdl.set('SpaceStep', 25);
mdl = mdl.set('Density',   1200);
mdl = mdl.set('Radius',    0.02);
mdl = mdl.set('Gravity',   [0,0,-9.81]);
mdl = mdl.set('E',         35);
mdl = mdl.set('Mu',        0.01);
mdl = mdl.set('Gain',      [K,0]);
mdl = mdl.set('Lambda',    K/X);
mdl = mdl.set('xia0',    [0,0,0,1,0,0].');

mdl = mdl.set('ActuationSpace',-1);

mdl = mdl.set('Point',...
    L0*[1,0,0,0,0.11,0.00,-0.01]);
end
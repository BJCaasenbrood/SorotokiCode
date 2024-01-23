clr;
%% Example: inverse kinematics solver for soft manipulator
% set parameters
[N,M] = deal(54,2);

% build polyspace
Y = chebyspace(N,M);

% construct Shapes class
shp = Shapes(Y,[0,M,M,0,0,0],'Quality',175);
shp = shp.setRamp(0.8);  % set ramp

[shp, q, g] = shp.solveIK([20,20,10]);
shp = shp.show();

view(30,30);
axis equal; axis tight;
box on; 
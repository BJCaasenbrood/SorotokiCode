clr; 
cd ..
cd ..
% assign parameters and build tensor
obj = preset.assets.tennis_racket;

[I1,I2,I3,m] = deal(18.36,40.32,44.517,0.1);
Mtt = blkdiag(diag([I1,I2,I3]),m * eye(3));

rgb = RigidBody(Mtt);
rgb.Gmodel = obj;

rgb.solver.TimeHorizon = 10;
rgb.solver.TimeStep = 1/30;

rgb.solver.sol.x(8)  = 1e-6;
rgb.solver.sol.x(9)  = 25;
rgb.solver.sol.x(12) = 2;

axis([-3,3,-3,30,-3,3]);
box on;
axis on;

while rgb.solver.Time < rgb.solver.TimeHorizon
    rgb = updateStatesFlow(rgb);
    rgb = rgb.render();

    drawnow limitrate;
    pause(1/60);
    view(30,30);
end

tcprintf('green', 'done! \n');

t = rgb.solver.sol.tout;
y = rgb.solver.sol.yout;

fig;
subplot(2,1,1);
plot(t,y(:,8),t,y(:,9),t,y(:,10));

subplot(2,1,2);
plot(t,y(:,8),t,y(:,9),t,y(:,10));

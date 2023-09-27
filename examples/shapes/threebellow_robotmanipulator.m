clr;
base = preset.shapes.threebellow;
[con1, con2] = preset.shapes.threebellow_connector;
fing1 = preset.shapes.threebellow_finger('phi',0);
fing2 = preset.shapes.threebellow_finger('phi',1);
fing3 = preset.shapes.threebellow_finger('phi',2);

base = base.addChild(fing1,1);
base = base.addChild(fing2,1);
base = base.addChild(fing3,1);
base = base.addChild(con1,0);
base = base.addChild(con2,1);

base.solver.sol.x(1) = 2.5e-2;
dt = 1/120;

base.show();
view(30,20);
axis tight;
drawnow;

while true

    % get time
    t = base.solver.Time;

    % update system and render
    % tic;
    base = base.update(dt);
    base = base.show();
    % dt = toc;
    drawnow limitrate;
    
end
